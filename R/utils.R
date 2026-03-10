# Level 1 function definitions -------------------------------------------------


#' Validate a GO graph object
#'
#' @description
#' Validates that an object is a \code{GO} or \code{GOSubgraph} instance
#' with the structural components required by GOcontext functions.
#'
#' @param go Object to validate.
#' @param arg \code{character(1)} Argument name for error reporting.
#' @param require_map \code{logical(1)} If \code{TRUE}, require an attached
#'   organism mapping in slot \code{map}.
#'
#' @return Invisibly \code{TRUE}.
#' @noRd
.assert_go_like <- function(
        go,
        arg = "go",
        require_map = FALSE
) {
    if (!methods::is(go, "GO") &&
        !methods::is(go, "GOSubgraph")) {
        rlang::abort(
            message = sprintf(
                "`%s` must be a GO or GOSubgraph object.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.character(go@ontology) ||
        length(go@ontology) != 1L ||
        is.na(go@ontology) ||
        !go@ontology %in% c("BP", "MF", "CC")) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `ontology`.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.character(go@version) ||
        length(go@version) != 1L ||
        is.na(go@version)) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `version`.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.data.frame(go@terms) ||
        !all(c("go_id", "term") %in% colnames(go@terms))) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `terms`.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.data.frame(go@edges) ||
        !all(c("child", "parent") %in% colnames(go@edges))) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `edges`.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.list(go@parents) || is.null(names(go@parents))) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `parents`.",
                arg
            ),
            arg = arg
        )
    }

    if (!is.list(go@children) || is.null(names(go@children))) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `children`.",
                arg
            ),
            arg = arg
        )
    }

    if (isTRUE(require_map)) {
        if (!is.data.frame(go@map) ||
            !all(c("go_id", "gene_id") %in% colnames(go@map)) ||
            nrow(go@map) == 0L) {
            rlang::abort(
                message = sprintf(
                    "`%s` must have a non-empty attached organism mapping.",
                    arg
                ),
                arg = arg
            )
        }
    }

    invisible(TRUE)
}


#' Induce an edge table on retained GO terms
#'
#' @description
#' Restricts an edge table to edges whose child and parent GO IDs are both
#' contained in the retained set of GO terms.
#'
#' @param edges \code{data.frame} Edge table with columns \code{child} and
#'   \code{parent}.
#' @param keep_ids \code{character()} vector of GO IDs to retain.
#'
#' @return A \code{data.frame} containing the induced edge set.
#'
#' @noRd
.induce_edges <- function(
        edges,
        keep_ids
) {
    if (!nrow(edges) || !length(keep_ids)) {
        return(edges[0, , drop = FALSE])
    }

    ok <- edges$child %in% keep_ids & edges$parent %in% keep_ids
    edges[ok, , drop = FALSE]
}


#' Induce adjacency lists on retained GO terms
#'
#' @description
#' Restricts the parent and child adjacency lists of a GO graph to a given
#' set of retained GO IDs.
#'
#' Adjacency entries pointing outside the retained node set are removed.
#'
#' @param parents Named list of parent relationships from the GO graph.
#' @param children Named list of child relationships from the GO graph.
#' @param keep_ids \code{character()} vector of GO IDs to retain.
#'
#' @return A named list with elements \code{parents} and \code{children},
#'   each restricted to \code{keep_ids}.
#'
#' @noRd
.induce_adjacency <- function(
        parents,
        children,
        keep_ids
) {
    parents_sub <- parents[keep_ids]
    children_sub <- children[keep_ids]

    parents_sub <- lapply(parents_sub, function(x) {
        if (is.null(x) || !length(x)) {
            return(character(0))
        }
        x <- x[x %in% keep_ids]
        if (length(x)) unique(x) else character(0)
    })

    children_sub <- lapply(children_sub, function(x) {
        if (is.null(x) || !length(x)) {
            return(character(0))
        }
        x <- x[x %in% keep_ids]
        if (length(x)) unique(x) else character(0)
    })

    names(parents_sub) <- keep_ids
    names(children_sub) <- keep_ids

    list(
        parents  = parents_sub,
        children = children_sub
    )
}


#' Identify GO terms passing gene set size thresholds
#'
#' @description
#' Determines which GO terms in a \code{GO} or \code{GOSubgraph} object
#' satisfy optional gene set size constraints based on the attached
#' organism mapping.
#'
#' The function first restricts the GO-to-gene mapping to terms present
#' in the current graph and then counts the number of mapped genes per
#' GO term.
#'
#' If \code{minGSSize} is provided, only terms with at least that many
#' mapped genes are retained. If \code{maxGSSize} is provided, only terms
#' with at most that many mapped genes are retained. If both are
#' provided, terms within the interval
#' \code{[minGSSize, maxGSSize]} are retained.
#'
#' This helper is used internally by \code{as_term2gene()} and
#' \code{as_term2name()} to ensure both exported mappings apply
#' identical gene set size filtering.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object with an attached
#'   organism mapping.
#' @param minGSSize \code{integer(1)} Minimum number of genes per GO term,
#'   or \code{NULL}.
#' @param maxGSSize \code{integer(1)} Maximum number of genes per GO term,
#'   or \code{NULL}.
#'
#' @return \code{character()} vector of GO IDs passing the size filter.
#'
#' @noRd
.filter_terms_by_size <- function(
        go,
        minGSSize,
        maxGSSize
) {
    map <- .restrict_go_map(
        map = go@map,
        go  = go
    )

    if (!nrow(map)) {
        return(character(0))
    }

    term_sizes <- table(map$go_id)
    keep <- rep(TRUE, length(term_sizes))

    if (!is.null(minGSSize)) {
        keep <- keep & term_sizes >= minGSSize
    }

    if (!is.null(maxGSSize)) {
        keep <- keep & term_sizes <= maxGSSize
    }

    names(term_sizes)[keep]
}


#' Validate gene set size thresholds
#'
#' @description
#' Validates the minimum and maximum gene set size parameters used for
#' optional filtering of GO terms by the number of mapped genes.
#'
#' If provided, \code{minGSSize} must be a positive scalar. If provided,
#' \code{maxGSSize} must be a positive scalar. If both are provided,
#' \code{maxGSSize} must be greater than or equal to \code{minGSSize}.
#'
#' @param minGSSize Minimum number of genes per GO term, or \code{NULL}.
#' @param maxGSSize Maximum number of genes per GO term, or \code{NULL}.
#'
#' @return Invisibly \code{TRUE} if validation succeeds.
#'
#' @noRd
.validate_gs_size <- function(
        minGSSize,
        maxGSSize
) {
    if (!is.null(minGSSize)) {
        if (!is.numeric(minGSSize)
            || length(minGSSize) != 1L
            || is.na(minGSSize)
            || minGSSize < 1) {
            rlang::abort(
                "`minGSSize` must be a positive integer or NULL.",
                arg = "minGSSize"
            )
        }
    }

    if (!is.null(maxGSSize)) {
        if (!is.numeric(maxGSSize)
            || length(maxGSSize) != 1L
            || is.na(maxGSSize)
            || maxGSSize < 1) {
            rlang::abort(
                "`maxGSSize` must be a positive integer or NULL.",
                arg = "maxGSSize"
            )
        }
    }

    if (!is.null(minGSSize) &&
        !is.null(maxGSSize) &&
        maxGSSize < minGSSize) {
        rlang::abort(
            "`maxGSSize` must be >= `minGSSize`.",
            arg = "maxGSSize"
        )
    }

    invisible(TRUE)
}
