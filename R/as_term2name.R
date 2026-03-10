#' Export GO term names as TERM2NAME
#'
#' @description
#' Returns a two-column \code{data.frame} suitable for
#' \code{clusterProfiler::enricher(TERM2NAME = ...)}.
#'
#' The returned table is restricted to GO terms present in the current
#' \code{GO} or \code{GOSubgraph} object.
#'
#' Optional size filtering can be applied based on the number of mapped
#' genes per GO term.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param minGSSize Minimum number of genes per GO term.
#' @param maxGSSize Maximum number of genes per GO term.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{name}.
#'
#' @examples
#' data("go_cc_ecoli", package = "GOcontext")
#'
#' # Export GO term names for the current graph
#' term2name <- as_term2name(go_cc_ecoli)
#' head(term2name)
#'
#' # Apply size filtering based on mapped genes
#' term2name_filtered <- as_term2name(
#'     go = go_cc_ecoli,
#'     minGSSize = 5,
#'     maxGSSize = 100
#' )
#' head(term2name_filtered)
#'
#' @export
as_term2name <- function(
        go,
        minGSSize = NULL,
        maxGSSize = NULL
) {
    use_size_filter <- !is.null(minGSSize) || !is.null(maxGSSize)
    .assert_go_like(
        go          = go,
        require_map = use_size_filter
    )

    terms <- .extract_term_metadata(go = go)

    if (!use_size_filter) {
        return(data.frame(
            term             = terms$go_id,
            name             = terms$term,
            stringsAsFactors = FALSE
        ))
    }

    .validate_gs_size(
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
    )
    keep_terms <- .filter_terms_by_size(
        go        = go,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
    )
    terms <- terms[terms$go_id %in% keep_terms, , drop = FALSE]

    data.frame(
        term             = terms$go_id,
        name             = terms$term,
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
    rownames(df) <- NULL
    df
}
