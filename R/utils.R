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

    if (!is.integer(go@depth)) {
        rlang::abort(
            message = sprintf(
                "`%s` has invalid slot `depth`.",
                arg
            ),
            arg = arg
        )
    }

    if (isTRUE(require_map)) {
        if (!is.data.frame(go@map) ||
            !all(c("go_id", "gene_id") %in% colnames(go@map))) {
            rlang::abort(
                message = sprintf(
                    "`%s` has no valid attached organism mapping. ",
                    arg
                ),
                arg = arg
            )
        }
    }

    invisible(TRUE)
}
