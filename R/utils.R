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
