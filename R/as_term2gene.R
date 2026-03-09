#' Export GO mappings as TERM2GENE
#'
#' @description
#' Returns a two-column \code{data.frame} suitable for
#' \code{clusterProfiler::enricher(TERM2GENE = ...)}.
#'
#' Requires an organism mapping attached via \code{attach_org()}.
#' The returned table is restricted to GO terms present in the current
#' \code{GO} or \code{GOSubgraph} object.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object with an attached
#'   organism mapping.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{gene}.
#'
#' @export
as_term2gene <- function(go) {
    .assert_go_like(
        go          = go,
        require_map = TRUE
    )
    map <- .restrict_go_map(
        map = go@map,
        go  = go
    )
    data.frame(
        term             = map$go_id,
        gene             = map$gene_id,
        stringsAsFactors = FALSE
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Restrict GO-to-gene mapping to terms present in the graph
#'
#' @description
#' Filters a GO-to-gene mapping table so that only mappings corresponding
#' to GO terms present in the supplied \code{GO} or \code{GOSubgraph}
#' object are retained.
#'
#' Rows with missing identifiers are removed and duplicate mappings are
#' dropped.
#'
#' @param map A \code{data.frame} containing GO-to-gene mappings with
#'   columns \code{go_id} and \code{gene_id}.
#' @param go A \code{GO} or \code{GOSubgraph} object used to define the
#'   set of valid GO terms.
#'
#' @return A filtered \code{data.frame} containing GO-to-gene mappings.
#'
#' @noRd
.restrict_go_map <- function(
        map,
        go
) {
    ids <- go@terms$go_id
    ids <- ids[!is.na(ids)]
    ids <- unique(ids)

    map <- map[!is.na(map$go_id) & !is.na(map$gene_id), , drop = FALSE]
    map <- map[map$go_id %in% ids, , drop = FALSE]
    map <- map[!duplicated(map), , drop = FALSE]
    rownames(map) <- NULL

    map
}
