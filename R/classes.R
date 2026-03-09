#' GO: full Gene Ontology DAG representation
#'
#' @description
#' Internal S4 class representing a full Gene Ontology (GO) directed
#' acyclic graph (DAG) for a single ontology (\code{"BP"}, \code{"MF"}, or
#' \code{"CC"}).
#'
#' The object stores term metadata and adjacency structures for efficient
#' graph traversal, subgraph induction, and organism-specific filtering.
#'
#' @slot ontology \code{character(1)} Ontology code (\code{"BP"},
#'   \code{"MF"}, or \code{"CC"}).
#' @slot version \code{character(1)} GO.db version string used to construct
#'   the graph.
#' @slot ids \code{character()} GO IDs present in the graph.
#' @slot terms \code{data.frame} Term-level metadata. Expected columns
#'   include at least \code{go_id}, \code{term}, \code{ontology}, and
#'   \code{obsolete}.
#' @slot edges \code{data.frame} Edge list with columns \code{child} and
#'   \code{parent}.
#' @slot parents \code{list} Named adjacency list: GO ID ->
#'   \code{character()} vector of direct parent GO IDs.
#' @slot children \code{list} Named adjacency list: GO ID ->
#'   \code{character()} vector of direct child GO IDs.
#' @slot map \code{data.frame} Attached organism mapping with columns
#'   \code{go_id} and \code{gene_id}. May be empty.
#'
#' @noRd
methods::setClass(
    "GO",
    slots = c(
        ontology = "character",
        version  = "character",
        ids      = "character",
        terms    = "data.frame",
        edges    = "data.frame",
        parents  = "list",
        children = "list",
        map      = "data.frame"
    )
)

setValidity("GO", function(object) {
    if (length(object@ontology) != 1L) {
        return("slot 'ontology' must be a scalar character.")
    }

    if (!object@ontology %in% c("BP", "MF", "CC")) {
        return("slot 'ontology' must be one of 'BP', 'MF', 'CC'.")
    }

    if (length(object@version) != 1L) {
        return("slot 'version' must be a scalar character.")
    }

    if (!identical(object@ids, object@terms$go_id)) {
        return("slot 'ids' must match terms$go_id.")
    }

    if (!is.data.frame(object@map)) {
        return("slot 'map' must be a data.frame.")
    }

    if (!all(c("go_id", "gene_id") %in% colnames(object@map))) {
        return("slot 'map' must contain columns 'go_id' and 'gene_id'.")
    }

    if (!identical(names(object@parents), object@ids)) {
        return("names(parents) must match ids.")
    }
    if (!identical(names(object@children), object@ids)) {
        return("names(children) must match ids.")
    }

    if (!all(c("child", "parent") %in% colnames(object@edges))) {
        return("slot 'edges' must contain columns 'child' and 'parent'.")
    }

    TRUE
})


#' GOSubgraph: restricted GO ontology subgraph
#'
#' @description
#' Internal S4 class representing a GO subgraph derived from a parent
#' \code{GO} or \code{GOSubgraph} object using
#' \code{subset_go(go, ids, mode = c("keep", "exclude"))}.
#'
#' @slot ontology \code{character(1)} Ontology code (\code{"BP"},
#'   \code{"MF"}, or \code{"CC"}).
#' @slot version \code{character(1)} GO.db version used to construct the
#'   graph.
#' @slot keep_ids \code{character()} GO IDs retained in the subgraph.
#' @slot drop_ids \code{character()} GO IDs removed from the parent graph.
#' @slot seed_ids \code{character()} GO IDs supplied to
#'   \code{subset_go()}.
#' @slot mode \code{character(1)} Subsetting mode, either \code{"keep"} or
#'   \code{"exclude"}.
#' @slot terms \code{data.frame} Term metadata for \code{keep_ids}.
#' @slot edges \code{data.frame} Induced edges with columns \code{child} and
#'   \code{parent}.
#' @slot parents \code{list} Adjacency: parents per node.
#' @slot children \code{list} Adjacency: children per node.
#' @slot map \code{data.frame} Attached organism mapping with columns
#'   \code{go_id} and \code{gene_id}. May be empty.
#'
#' @noRd
methods::setClass(
    "GOSubgraph",
    slots = c(
        ontology = "character",
        version  = "character",
        keep_ids = "character",
        drop_ids = "character",
        seed_ids = "character",
        mode     = "character",
        terms    = "data.frame",
        edges    = "data.frame",
        parents  = "list",
        children = "list",
        map      = "data.frame"
    )
)

setValidity("GOSubgraph", function(object) {
    if (length(object@ontology) != 1L) {
        return("slot 'ontology' must be a scalar character.")
    }

    if (!object@ontology %in% c("BP", "MF", "CC")) {
        return("slot 'ontology' must be one of 'BP', 'MF', 'CC'.")
    }

    if (length(object@version) != 1L) {
        return("slot 'version' must be a scalar character.")
    }

    if (length(object@mode) != 1L) {
        return("slot 'mode' must be a scalar character.")
    }

    if (!object@mode %in% c("keep", "exclude")) {
        return("slot 'mode' must be one of 'keep' or 'exclude'.")
    }

    if (!identical(object@keep_ids, object@terms$go_id)) {
        return("slot 'keep_ids' must match terms$go_id.")
    }

    if (!is.data.frame(object@map)) {
        return("slot 'map' must be a data.frame.")
    }

    if (!all(c("go_id", "gene_id") %in% colnames(object@map))) {
        return("slot 'map' must contain columns 'go_id' and 'gene_id'.")
    }

    if (!identical(names(object@parents), object@keep_ids)) {
        return("names(parents) must match keep_ids.")
    }
    if (!identical(names(object@children), object@keep_ids)) {
        return("names(children) must match keep_ids.")
    }

    if (!all(c("child", "parent") %in% colnames(object@edges))) {
        return("slot 'edges' must contain columns 'child' and 'parent'.")
    }

    TRUE
})
