#' Attach organism GO mappings to a GO graph
#'
#' @description
#' Retrieves GO-to-gene mappings from an \code{OrgDb} and attaches them to a
#' \code{GO} or \code{GOSubgraph} object for downstream filtering and export.
#'
#' The mapping is restricted to GO IDs present in \code{go}. No graph
#' restriction is performed here; use \code{filter_mapped()} for that.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param OrgDb An \code{OrgDb} object (e.g. \code{org.EcK12.eg.db}).
#' @param keytype \code{character(1)} Gene identifier type to use from the
#'   \code{OrgDb}. Default is \code{"ENTREZID"}.
#'
#' @return The same object as \code{go}, with an attached mapping table.
#'
#' @examples
#' data("go_cc", package = "GOcontext")
#'
#' if (requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
#'     go_cc_ecoli <- attach_org(
#'         go = go_cc,
#'         OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
#'         keytype = "ENTREZID"
#'     )
#'
#'     go_cc_ecoli
#' }
#'
#' @export
attach_org <- function(
        go,
        OrgDb,
        keytype = "ENTREZID"
) {
    keytype <- trimws(keytype)
    .attach_org_validate_inputs(
        go      = go,
        OrgDb   = OrgDb,
        keytype = keytype
    )

    map <- .attach_org_build_map(
        go      = go,
        OrgDb   = OrgDb,
        keytype = keytype
    )

    .attach_org_set_map(
        go  = go,
        map = map
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Validate inputs for organism attachment
#'
#' @description
#' Validates the inputs provided to \code{attach_org()}. Ensures that the GO
#' object is a valid \code{GO} or \code{GOSubgraph}, that the supplied
#' \code{OrgDb} is a valid annotation database, and that the requested
#' \code{keytype} exists in the database.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param OrgDb An \code{OrgDb} object providing organism annotations.
#' @param keytype \code{character(1)} Gene identifier type used for mapping.
#'
#' @return Invisibly \code{TRUE} if validation succeeds.
#'
#' @noRd
.attach_org_validate_inputs <- function(
        go,
        OrgDb,
        keytype
) {
    .assert_go_like(
        go  = go,
        arg = "go"
    )

    if (!methods::is(OrgDb, "OrgDb")) {
        rlang::abort(
            message = "`OrgDb` must be an OrgDb object.",
            arg = "OrgDb"
        )
    }

    if (!is.character(keytype) ||
        length(keytype) != 1L ||
        is.na(keytype)) {
        rlang::abort(
            message = "`keytype` must be a non-NA character(1).",
            arg = "keytype"
        )
    }

    if (!nzchar(keytype)) {
        rlang::abort(
            message = "`keytype` must not be empty.",
            arg = "keytype"
        )
    }

    available_keytypes <- .annot_keytypes(OrgDb)
    if (!keytype %in% available_keytypes) {
        rlang::abort(
            message = sprintf(
                paste0(
                    "`keytype` (%s) is not available in `OrgDb`. ",
                    "Available keytypes: %s"
                ),
                sQuote(keytype),
                paste(available_keytypes, collapse = ", ")
            ),
            arg = "keytype"
        )
    }

    available_columns <- .annot_columns(OrgDb)
    if (!keytype %in% available_columns) {
        rlang::abort(
            message = sprintf(
                paste(
                    "`keytype` (%s) is not available as a selectable",
                    "column in `OrgDb`."
                ),
                sQuote(keytype)
            ),
            arg = "keytype"
        )
    }

    invisible(TRUE)
}


#' Build GO-to-gene mapping from an OrgDb
#'
#' @description
#' Retrieves GO-to-gene mappings from an \code{OrgDb} and restricts them to
#' GO terms present in the supplied GO graph.
#'
#' The function prefers a GO-driven lookup, using the GO IDs already present
#' in \code{go} as query keys. If the supplied \code{OrgDb} does not support
#' \code{"GO"} as a keytype, it falls back to retrieving all keys for the
#' requested gene identifier type and filtering the result to the GO graph.
#'
#' The resulting mapping contains two columns: \code{go_id} and
#' \code{gene_id}.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param OrgDb An \code{OrgDb} object providing organism annotations.
#' @param keytype \code{character(1)} Gene identifier type used for mapping.
#'
#' @return A \code{data.frame} with columns \code{go_id} and \code{gene_id}.
#'
#' @noRd
.attach_org_build_map <- function(
        go,
        OrgDb,
        keytype
) {
    go_ids <- unique(go@terms$go_id)
    go_ids <- go_ids[!is.na(go_ids)]

    if (!length(go_ids)) {
        return(.empty_go_map())
    }

    available_keytypes <- .annot_keytypes(OrgDb)
    available_columns  <- .annot_columns(OrgDb)

    use_go_keytype <-
        "GO" %in% available_keytypes && keytype %in% available_columns

    if (use_go_keytype) {
        res <- .silent_annotationdbi_select(
            x       = OrgDb,
            keys    = go_ids,
            keytype = "GO",
            columns = keytype
        )
    } else {
        keys <- .annot_keys(
            OrgDb = OrgDb,
            keytype = keytype
        )

        if (!length(keys)) {
            rlang::abort(
                message = sprintf(
                    "No keys returned from `OrgDb` for keytype %s.",
                    sQuote(keytype)
                )
            )
        }

        res <- .silent_annotationdbi_select(
            x       = OrgDb,
            keys    = keys,
            keytype = keytype,
            columns = "GO"
        )
    }

    if (!nrow(res)) {
        return(.empty_go_map())
    }

    if (!all(c("GO", keytype) %in% colnames(res))) {
        rlang::abort(
            message = paste(
                "`AnnotationDbi::select()` did not return the expected",
                "columns for GO mapping."
            )
        )
    }

    res <- res[!is.na(res$GO), , drop = FALSE]
    res <- res[!is.na(res[[keytype]]), , drop = FALSE]
    res <- res[res$GO %in% go_ids, , drop = FALSE]

    if (!nrow(res)) {
        return(.empty_go_map())
    }

    map <- data.frame(
        go_id            = as.character(res$GO),
        gene_id          = as.character(res[[keytype]]),
        stringsAsFactors = FALSE
    )

    map <- map[!duplicated(map), , drop = FALSE]
    rownames(map) <- NULL

    map
}


#' Attach a GO-to-gene mapping to a GO object
#'
#' @description
#' Stores the GO-to-gene mapping on the supplied GO object by writing the
#' mapping table into the \code{map} slot.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param map A \code{data.frame} containing GO-to-gene mappings.
#'
#' @return The same GO object with the mapping attached.
#'
#' @noRd
.attach_org_set_map <- function(
        go,
        map
) {
    methods::slot(go, "map") <- map
    go
}


# Level 2 helpers --------------------------------------------------------------


#' Construct an empty GO-to-gene mapping table
#'
#' @return A zero-row \code{data.frame} with columns \code{go_id} and
#'   \code{gene_id}.
#' @noRd
.empty_go_map <- function() {
    data.frame(
        go_id            = character(0),
        gene_id          = character(0),
        stringsAsFactors = FALSE
    )
}


.annot_keytypes <- function(OrgDb) {
    AnnotationDbi::keytypes(OrgDb)
}


.annot_columns <- function(OrgDb) {
    AnnotationDbi::columns(OrgDb)
}


.annot_keys <- function(
        OrgDb,
        keytype
) {
    AnnotationDbi::keys(
        x = OrgDb,
        keytype = keytype
    )
}
