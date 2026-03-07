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
#' @export
filter_mapped <- function(go) {
    .assert_go_like(
        go = go,
        require_map = TRUE
    )

    map <- go@map
    universe <- .extract_universe(go = go)
    keep_ids <- .find_keep_ids(
        map = map,
        universe = universe
        )

    edges_sub <- .induce_edges(
        edges = go@edges,
        keep_ids = keep_ids
    )

    adj <- .induce_adj(
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
        seed_ids = character(0),
        mode     = "keep",
        terms    = terms_sub,
        edges    = edges_sub,
        parents  = adj$parents,
        children = adj$children,
        depth    = integer(0),
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
    intersect(
        keep,
        universe
        )
}


#' Induce an edge table on mapped GO terms
#'
#' @description
#' Restricts an edge table to edges whose child and parent GO IDs are both
#' contained in the retained set of mapped GO terms.
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
        return(data.frame(
            child = character(0),
            parent = character(0),
            stringsAsFactors = FALSE
        ))
    }

    ok <- edges$child %in% keep_ids & edges$parent %in% keep_ids
    edges[ok, , drop = FALSE]
}


#' Induce adjacency lists on mapped GO terms
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
.induce_adj <- function(
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
