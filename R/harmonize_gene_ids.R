#' Harmonize gene identifiers using an OrgDb annotation package
#'
#' @description
#' Harmonizes a character vector of gene identifiers to a target identifier
#' type using an \code{OrgDb} object.
#'
#' The function is intentionally strict and follows a "don't guess"
#' strategy. Input identifiers are kept unchanged when they already match
#' \code{to_keytype}. Otherwise, candidate source identifier types in
#' \code{from_keytypes} are checked.
#'
#' Harmonization is only performed when exactly one source identifier type
#' matches and this yields exactly one target identifier. All non-unique or
#' unresolved cases are returned as \code{NA} in the \code{harmonized}
#' column.
#'
#' Missing input values are accepted and returned unchanged in the output
#' structure.
#'
#' The returned \code{data.frame} is always row-aligned with \code{genes},
#' including duplicated and missing values.
#'
#' @param genes A character vector of gene identifiers to harmonize.
#' @param OrgDb An \code{OrgDb} object providing organism-specific gene
#'   identifier mappings.
#' @param to_keytype A non-empty character scalar giving the target
#'   identifier type.
#' @param from_keytypes A character vector of candidate source identifier
#'   types to try when an input value does not already match
#'   \code{to_keytype}. Duplicates are removed internally, and
#'   \code{to_keytype} is excluded if present.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{\code{input}}{The original input identifier.}
#'   \item{\code{matches_target}}{Logical indicator of whether the input
#'   already matched \code{to_keytype}.}
#'   \item{\code{from_match_status}}{Status of matching against
#'   \code{from_keytypes}. One of \code{"none"}, \code{"unique"},
#'   \code{"ambiguous"}, or \code{NA} when the input already matched
#'   \code{to_keytype}.}
#'   \item{\code{to_match_status}}{Status of mapping from the uniquely
#'   matched source identifier type to \code{to_keytype}. One of
#'   \code{"none"}, \code{"unique"}, \code{"ambiguous"}, or \code{NA}
#'   when no target mapping was attempted.}
#'   \item{\code{harmonized}}{The harmonized identifier, or \code{NA} when
#'   no unique harmonization was possible.}
#' }
#'
#' @examples
#' if (requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
#'     genes <- c(
#'         "b0002",    # SYMBOL
#'         "thrA",     # ALIAS
#'         "944742",   # ENTREZID
#'         "not_a_gene",
#'         NA_character_
#'     )
#'
#'     res <- harmonize_gene_ids(
#'         genes = genes,
#'         OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
#'         to_keytype = "ENTREZID",
#'         from_keytypes = c("SYMBOL", "ALIAS")
#'     )
#'
#'     head(res)
#'
#'     # Extract harmonized identifiers
#'     res$harmonized
#' }
#'
#' @export
harmonize_gene_ids <- function(
        genes,
        OrgDb,
        to_keytype,
        from_keytypes = c("ALIAS", "ENTREZID", "ENSEMBL")
) {
    .harmonize_gene_ids_validate_inputs(
        genes = genes,
        OrgDb = OrgDb,
        to_keytype = to_keytype,
        from_keytypes = from_keytypes
    )

    from_keytypes <- .normalize_from_keytypes(
        from_keytypes = from_keytypes,
        to_keytype = to_keytype
    )

    is_missing <- is.na(genes)
    genes_work <- genes[!is_missing]

    key_sets <- .build_key_sets(
        OrgDb = OrgDb,
        to_keytype = to_keytype,
        from_keytypes = from_keytypes
    )

    direct_match <- genes_work %in% key_sets$target_keys

    from_info <- .classify_from_matches(
        genes = genes_work,
        direct_match = direct_match,
        from_keytypes = from_keytypes,
        from_key_sets = key_sets$from_key_sets
    )

    to_info <- .map_to_target_batch(
        genes = genes_work,
        matched_from_keytype = from_info$matched_from_keytype,
        OrgDb = OrgDb,
        to_keytype = to_keytype
    )

    results <- data.frame(
        input = genes,
        matches_target = FALSE,
        from_match_status = NA_character_,
        to_match_status = NA_character_,
        harmonized = NA_character_,
        stringsAsFactors = FALSE
    )

    idx <- which(!is_missing)
    idx_direct <- idx[direct_match]

    results$matches_target[idx] <- direct_match
    results$from_match_status[idx] <- from_info$from_match_status
    results$to_match_status[idx] <- to_info$to_match_status
    results$harmonized[idx] <- to_info$harmonized

    results$from_match_status[idx_direct] <- NA_character_
    results$to_match_status[idx_direct] <- NA_character_
    results$harmonized[idx_direct] <- genes[idx_direct]

    results
}


# Level 1 helpers --------------------------------------------------------------


#' Validate inputs for gene identifier harmonization
#'
#' @description
#' Performs input validation for \code{harmonize_gene_ids()}.
#'
#' The function checks that the supplied arguments have the expected
#' types and values and that the requested identifier types exist in the
#' provided \code{OrgDb} object.
#'
#' Specifically, it verifies that:
#' \itemize{
#'   \item \code{OrgDb} is an \code{OrgDb} object.
#'   \item \code{genes} is a character vector.
#'   \item \code{to_keytype} is a non-empty character scalar.
#'   \item \code{from_keytypes}, if provided, is a character vector without
#'   \code{NA} or empty strings.
#'   \item \code{to_keytype} and all \code{from_keytypes} exist among the
#'   available keytypes in \code{OrgDb}.
#' }
#'
#' Errors are raised using \code{rlang::abort()} when validation fails.
#'
#' @param genes A character vector of gene identifiers.
#' @param OrgDb An \code{OrgDb} object providing organism-specific
#'   identifier mappings.
#' @param to_keytype A character scalar specifying the target identifier
#'   type.
#' @param from_keytypes A character vector of candidate source identifier
#'   types or \code{NULL}.
#'
#' @return Invisibly returns \code{NULL}. The function is called only for
#' argument validation and produces an error if checks fail.
#'
#' @noRd
.harmonize_gene_ids_validate_inputs <- function(
        genes,
        OrgDb,
        to_keytype,
        from_keytypes
) {
    if (!methods::is(OrgDb, "OrgDb")) {
        rlang::abort(
            message = "`OrgDb` must be an OrgDb object.",
            arg = "OrgDb"
        )
    }

    if (!is.character(genes)) {
        rlang::abort(
            message = "`genes` must be a character vector.",
            arg = "genes"
        )
    }

    if (!is.character(to_keytype) ||
        length(to_keytype) != 1L ||
        is.na(to_keytype) ||
        !nzchar(to_keytype)) {
        rlang::abort(
            message = "`to_keytype` must be a non-empty character(1).",
            arg = "to_keytype"
        )
    }

    if (!is.null(from_keytypes) &&
        (!is.character(from_keytypes) ||
         anyNA(from_keytypes) ||
         any(!nzchar(from_keytypes)))) {
        rlang::abort(
            message = paste(
                "`from_keytypes` must be NULL or a character vector",
                "containing no NA or empty strings."
            ),
            arg = "from_keytypes"
        )
    }

    available_keytypes <- AnnotationDbi::keytypes(OrgDb)

    if (!to_keytype %in% available_keytypes) {
        rlang::abort(
            message = sprintf(
                paste0(
                    "`to_keytype` (%s) is not available in `OrgDb`. ",
                    "Available keytypes: %s"
                ),
                sQuote(to_keytype),
                paste(available_keytypes, collapse = ", ")
            ),
            arg = "to_keytype"
        )
    }

    if (!is.null(from_keytypes)) {
        bad_from <- setdiff(unique(from_keytypes), available_keytypes)

        if (length(bad_from) > 0L) {
            rlang::abort(
                message = sprintf(
                    paste0(
                        "The following `from_keytypes` are not ",
                        "available in `OrgDb`: %s"
                    ),
                    paste(bad_from, collapse = ", ")
                ),
                arg = "from_keytypes"
            )
        }
    }

    invisible(NULL)
}


#' Normalize candidate source keytypes
#'
#' @description
#' Cleans and standardizes the \code{from_keytypes} argument used in
#' \code{harmonize_gene_ids()}.
#'
#' If \code{from_keytypes} is \code{NULL}, it is replaced with an empty
#' character vector. Duplicate entries are removed, and the target
#' identifier type \code{to_keytype} is excluded if present.
#'
#' This ensures that candidate source keytypes are unique and do not
#' include the target identifier type.
#'
#' @param from_keytypes A character vector of candidate source identifier
#'   types, or \code{NULL}.
#' @param to_keytype A character scalar giving the target identifier type.
#'
#' @return A character vector of normalized source keytypes.
#'
#' @noRd
.normalize_from_keytypes <- function(
        from_keytypes,
        to_keytype
) {
    if (is.null(from_keytypes)) {
        from_keytypes <- character(0)
    }

    from_keytypes <- unique(from_keytypes)
    from_keytypes <- setdiff(from_keytypes, to_keytype)

    from_keytypes
}


#' Build key lookup sets for identifier matching
#'
#' @description
#' Retrieves identifier sets from an \code{OrgDb} object used for
#' matching and harmonization in \code{harmonize_gene_ids()}.
#'
#' The function extracts all identifiers available for the target
#' \code{to_keytype} as well as identifier sets for each candidate
#' source type in \code{from_keytypes}. These sets are used for fast
#' membership testing during the harmonization process.
#'
#' @param OrgDb An \code{OrgDb} object providing organism-specific
#'   identifier mappings.
#' @param to_keytype A character scalar specifying the target identifier
#'   type.
#' @param from_keytypes A character vector of candidate source identifier
#'   types.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{target_keys}}{A character vector containing all
#'   identifiers available for \code{to_keytype}.}
#'   \item{\code{from_key_sets}}{A named list of character vectors,
#'   where each element contains all identifiers available for the
#'   corresponding source keytype.}
#' }
#'
#' @noRd
.build_key_sets <- function(
        OrgDb,
        to_keytype,
        from_keytypes
) {
    target_keys <- unique(as.character(AnnotationDbi::keys(
        x = OrgDb,
        keytype = to_keytype
    )))

    from_key_sets <- lapply(
        from_keytypes,
        function(kt) {
            unique(as.character(AnnotationDbi::keys(
                x = OrgDb,
                keytype = kt
            )))
        }
    )
    names(from_key_sets) <- from_keytypes

    list(
        target_keys = target_keys,
        from_key_sets = from_key_sets
    )
}


#' Classify matches against candidate source keytypes
#'
#' @description
#' Determines, for each input identifier, whether it matches none, one, or
#' multiple candidate source keytypes in \code{from_keytypes}.
#'
#' Direct matches to the target identifier type are excluded from this
#' classification step and retain \code{NA} in the returned
#' \code{from_match_status}.
#'
#' For inputs with exactly one source-side match, the matching source
#' keytype is also recorded.
#'
#' @param genes A character vector of input identifiers to classify.
#' @param direct_match A logical vector indicating which inputs already
#'   match the target identifier type.
#' @param from_keytypes A character vector of candidate source identifier
#'   types.
#' @param from_key_sets A named list of character vectors containing the
#'   available identifiers for each source keytype.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{from_match_status}}{A character vector with values
#'   \code{"none"}, \code{"unique"}, \code{"ambiguous"}, or \code{NA}
#'   for direct target matches.}
#'   \item{\code{matched_from_keytype}}{A character vector giving the
#'   uniquely matched source keytype, or \code{NA} otherwise.}
#' }
#'
#' @noRd
.classify_from_matches <- function(
        genes,
        direct_match,
        from_keytypes,
        from_key_sets
) {
    n_genes <- length(genes)

    from_match_status <- rep(NA_character_, n_genes)
    matched_from_keytype <- rep(NA_character_, n_genes)

    idx <- which(!direct_match)

    if (length(idx) == 0L || length(from_keytypes) == 0L) {
        from_match_status[idx] <- "none"

        return(list(
            from_match_status = from_match_status,
            matched_from_keytype = matched_from_keytype
        ))
    }

    hit_matrix <- vapply(
        from_key_sets,
        function(keys_i) {
            genes[idx] %in% keys_i
        },
        logical(length(idx))
    )

    if (is.null(dim(hit_matrix))) {
        hit_matrix <- matrix(
            hit_matrix,
            nrow = length(idx),
            ncol = length(from_key_sets)
        )
        colnames(hit_matrix) <- names(from_key_sets)
    }

    hit_counts <- rowSums(hit_matrix)

    from_match_status[idx][hit_counts == 0L] <- "none"
    from_match_status[idx][hit_counts == 1L] <- "unique"
    from_match_status[idx][hit_counts > 1L] <- "ambiguous"

    unique_idx <- idx[hit_counts == 1L]

    if (length(unique_idx) > 0L) {
        unique_cols <- max.col(hit_matrix[hit_counts == 1L, , drop = FALSE])
        matched_from_keytype[unique_idx] <- colnames(hit_matrix)[unique_cols]
    }

    list(
        from_match_status = from_match_status,
        matched_from_keytype = matched_from_keytype
    )
}


#' Map uniquely matched source identifiers to the target keytype
#'
#' @description
#' Performs batched mapping from source identifier types to the target
#' identifier type in \code{harmonize_gene_ids()}.
#'
#' Inputs are grouped by their uniquely matched source keytype, and each
#' group is queried in a single call to \code{AnnotationDbi::select()}.
#' This improves performance compared with querying each identifier
#' individually.
#'
#' For each input identifier with a unique source keytype, the function
#' determines whether mapping to \code{to_keytype} yields no match, a
#' unique match, or multiple matches. Only unique target matches are
#' returned in \code{harmonized}.
#'
#' @param genes A character vector of input identifiers.
#' @param matched_from_keytype A character vector giving the uniquely
#'   matched source keytype for each input identifier, or \code{NA} when
#'   no unique source keytype was identified.
#' @param OrgDb An \code{OrgDb} object providing organism-specific
#'   identifier mappings.
#' @param to_keytype A character scalar specifying the target identifier
#'   type.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{to_match_status}}{A character vector with values
#'   \code{"none"}, \code{"unique"}, \code{"ambiguous"}, or \code{NA}
#'   when no target-side mapping was attempted.}
#'   \item{\code{harmonized}}{A character vector of harmonized target
#'   identifiers, or \code{NA} when no unique target mapping was
#'   obtained.}
#' }
#'
#' @noRd
.map_to_target_batch <- function(
        genes,
        matched_from_keytype,
        OrgDb,
        to_keytype
) {
    n_genes <- length(genes)

    to_match_status <- rep(NA_character_, n_genes)
    harmonized <- rep(NA_character_, n_genes)

    idx <- which(!is.na(matched_from_keytype))

    if (length(idx) == 0L) {
        return(list(
            to_match_status = to_match_status,
            harmonized = harmonized
        ))
    }

    split_idx <- split(idx, matched_from_keytype[idx])

    for (kt in names(split_idx)) {
        idx_kt <- split_idx[[kt]]
        genes_kt <- genes[idx_kt]

        sel <- .silent_annotationdbi_select(
            x       = OrgDb,
            keys    = genes_kt,
            keytype = kt,
            columns = to_keytype
        )

        if (!nrow(sel) || !all(c(kt, to_keytype) %in% names(sel))) {
            to_match_status[idx_kt] <- "none"
            next
        }

        sel[[kt]] <- as.character(sel[[kt]])
        sel[[to_keytype]] <- as.character(sel[[to_keytype]])

        target_list <- split(sel[[to_keytype]], sel[[kt]])

        for (i in idx_kt) {
            vals <- target_list[[genes[[i]]]]

            if (is.null(vals)) {
                to_match_status[[i]] <- "none"
                next
            }

            vals <- unique(vals[!is.na(vals)])

            if (length(vals) == 0L) {
                to_match_status[[i]] <- "none"
                next
            }

            if (length(vals) > 1L) {
                to_match_status[[i]] <- "ambiguous"
                next
            }

            to_match_status[[i]] <- "unique"
            harmonized[[i]] <- vals
        }
    }

    list(
        to_match_status = to_match_status,
        harmonized = harmonized
    )
}
