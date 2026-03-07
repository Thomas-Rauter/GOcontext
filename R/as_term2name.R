#' Export GO term names as TERM2NAME
#'
#' @description
#' Returns a two-column \code{data.frame} suitable for
#' \code{clusterProfiler::enricher(TERM2NAME = ...)}.
#'
#' The returned table is restricted to GO terms present in the current
#' \code{GO} or \code{GOSubgraph} object.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{name}.
#'
#' @export
as_term2name <- function(go) {
    .assert_go_like(go = go)

    terms <- .extract_term_metadata(go = go)

    data.frame(
        term = terms$go_id,
        name = terms$term,
        stringsAsFactors = FALSE
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Extract GO term metadata from a GO object
#'
#' @description
#' Retrieves the GO term identifiers and term names from a \code{GO} or
#' \code{GOSubgraph} object and returns a clean two-column table used by
#' \code{as_term2name()}.
#'
#' Rows with missing identifiers or names are removed, duplicated GO IDs
#' are dropped, and the result is ordered by GO ID.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#'
#' @return A \code{data.frame} with columns \code{go_id} and \code{term}.
#'
#' @noRd
.extract_term_metadata <- function(go) {
    df <- go@terms[, c("go_id", "term"), drop = FALSE]
    df <- df[!is.na(df$go_id) & !is.na(df$term), , drop = FALSE]
    df <- df[!duplicated(df$go_id), , drop = FALSE]
    df <- df[order(df$go_id), , drop = FALSE]
    rownames(df) <- NULL

    df
}
