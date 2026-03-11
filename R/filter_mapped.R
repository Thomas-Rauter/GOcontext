#' Restrict a GO graph to mapped GO terms
#'
#' @description
#' Restricts a \code{GO} or \code{GOSubgraph} object to GO terms that have
#' at least one mapped gene in the attached organism annotation.
#'
#' The mapping must have been attached previously using
#' \code{attach_org()}. The mapping itself is not modified; only the
#' graph structure is restricted.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object with an attached
#'   organism mapping.
#'
#' @return A \code{GOSubgraph} restricted to GO terms with at least one
#'   mapped gene.
#'
#' @examples
#' data("go_cc_ecoli", package = "GOcontext")
#'
#' go_cc_ecoli_mapped <- filter_mapped(go_cc_ecoli)
#' go_cc_ecoli_mapped
#'
#' @export
filter_mapped <- function(go) {
    .assert_go_like(
        go = go,
        require_map = TRUE
    )

    map <- go@map
    universe <- .extract_universe(go = go)

    mapped_ids <- .find_keep_ids(
        map = map,
        universe = universe
    )

    keep_ids <- .ancestor_closure(
        seed_ids = mapped_ids,
        parents  = go@parents,
        universe = universe
    )

    edges_sub <- .induce_edges(
        edges = go@edges,
        keep_ids = keep_ids
    )

    adj <- .induce_adjacency(
        parents = go@parents,
        children = go@children,
        keep_ids = keep_ids
    )

    terms_sub <- .subset_to_mapped(
        go = go,
        keep_ids = keep_ids
    )

    drop_ids <- setdiff(
        universe,
        keep_ids
    )

    methods::new(
        "GOSubgraph",
        ontology = go@ontology,
        version  = go@version,
        keep_ids = keep_ids,
        drop_ids = drop_ids,
        seed_ids = mapped_ids,
        mode     = "keep",
        terms    = terms_sub,
        edges    = edges_sub,
        parents  = adj$parents,
        children = adj$children,
        map      = go@map
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Extract the GO ID universe from a GO graph
#'
#' @description
#' Retrieves the set of GO identifiers present in a \code{GO} or
#' \code{GOSubgraph} object from its term metadata table.
#'
#' Missing GO IDs are removed and the result is returned as a unique
#' character vector.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#'
#' @return \code{character()} vector of GO IDs present in the graph.
#'
#' @noRd
.extract_universe <- function(go) {
    ids <- go@terms$go_id
    ids <- ids[!is.na(ids)]
    unique(ids)
}


#' Identify GO terms with organism mappings
#'
#' @description
#' Determines which GO terms in the current graph have at least one mapped
#' gene in the attached organism annotation.
#'
#' The returned GO IDs are restricted to the supplied graph universe.
#'
#' @param map A \code{data.frame} containing GO-to-gene mappings with
#'   columns \code{go_id} and \code{gene_id}.
#' @param universe \code{character()} vector of GO IDs present in the
#'   current graph.
#'
#' @return \code{character()} vector of GO IDs with at least one mapping.
#'
#' @noRd
.find_keep_ids <- function(
        map,
        universe
) {
    keep <- unique(map$go_id)
    keep <- keep[!is.na(keep)]
    universe[universe %in% keep]
    intersect(
        keep,
        universe
        )
}


#' Compute ancestor closure of GO terms
#'
#' @description
#' Starting from a set of seed GO IDs, traverses the GO DAG upward through
#' parent relationships and returns the seed terms together with all of
#' their ancestors.
#'
#' The result is optionally restricted to a supplied universe of GO IDs.
#'
#' @param seed_ids `character()` GO IDs to start from.
#' @param parents Named `list` mapping each GO ID to its direct parent GO IDs.
#' @param universe Optional `character()` vector of GO IDs allowed in the
#'   output. If supplied, the returned IDs are restricted to this set.
#'
#' @return `character()` vector containing `seed_ids` and all of their
#'   ancestors.
#'
#' @noRd
.ancestor_closure <- function(
        seed_ids,
        parents,
        universe = NULL
) {
    seed_ids <- unique(seed_ids)
    seed_ids <- seed_ids[!is.na(seed_ids)]

    if (!is.null(universe)) {
        .validate_keep_ids(universe, arg = "universe")
        seed_ids <- intersect(seed_ids, universe)
    }

    if (!length(seed_ids)) {
        return(character(0))
    }

    keep <- seed_ids
    queue <- seed_ids

    while (length(queue)) {
        node <- queue[[1L]]
        queue <- queue[-1L]

        node_parents <- parents[[node]]

        if (is.null(node_parents)) {
            next
        }

        node_parents <- unique(node_parents)
        node_parents <- node_parents[!is.na(node_parents)]

        if (!is.null(universe)) {
            node_parents <- intersect(node_parents, universe)
        }

        new_parents <- setdiff(node_parents, keep)

        if (length(new_parents)) {
            keep <- c(keep, new_parents)
            queue <- c(queue, new_parents)
        }
    }

    unique(keep)
}


#' Subset GO term metadata to mapped GO terms
#'
#' @description
#' Extracts the term metadata corresponding to the retained GO IDs after
#' organism-based filtering.
#'
#' If no GO IDs are retained, an empty term table with the same columns is
#' returned.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param keep_ids \code{character()} vector of GO IDs to retain.
#'
#' @return A \code{data.frame} containing term metadata for
#'   \code{keep_ids}.
#'
#' @noRd
.subset_to_mapped <- function(
        go,
        keep_ids
) {
    if (!length(keep_ids)) {
        return(go@terms[0, , drop = FALSE])
    }

    terms_sub <- go@terms[match(keep_ids, go@terms$go_id), , drop = FALSE]
    terms_sub <- terms_sub[!is.na(terms_sub$go_id), , drop = FALSE]
    terms_sub
}
