#' Restrict a GO graph by keeping or excluding selected branches
#'
#' @description
#' Given one or more GO IDs, induce a subgraph by either keeping those terms
#' and all their descendants, or excluding those terms and all their
#' descendants, within the current graph.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param ids \code{character()} Vector of GO IDs used as subsetting seeds.
#' @param mode \code{character(1)} One of \code{"keep"} or
#'   \code{"exclude"}.
#'
#' @return A \code{GOSubgraph} object.
#'
#' @export
subset_go <- function(
        go,
        ids,
        mode = c("keep", "exclude")
) {
    .subset_go_validate_inputs(
        go = go,
        ids = ids
        )
    mode <- match.arg(mode)

    universe <- .extract_go_ids(go)
    ids <- unique(ids)

    .check_ids_in_graph(
        ids = ids,
        universe = universe
        )

    seed_desc <- .expand_to_descendants(
        seeds = ids,
        children = go@children,
        universe = universe
    )

    keep_ids <- .resolve_retained_ids(
        universe = universe,
        seed_desc = seed_desc,
        mode = mode
    )

    edges_sub <- .induce_edges(
        edges = go@edges,
        keep_ids = keep_ids
        )

    adj_sub <- .induce_adjacency(
        parents = go@parents,
        children = go@children,
        keep_ids = keep_ids
    )

    terms_sub <- .restrict_terms(
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
        seed_ids = ids,
        mode     = mode,
        terms    = terms_sub,
        edges    = edges_sub,
        parents  = adj_sub$parents,
        children = adj_sub$children,
        map      = go@map
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Validate inputs
#'
#' @description
#' Validates the GO graph object and the vector of seed GO IDs supplied for
#' graph restriction.
#'
#' The GO object must be a valid \code{GO} or \code{GOSubgraph}, and
#' \code{ids} must be a non-empty character vector without missing values.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param ids \code{character()} vector of seed GO IDs.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @noRd
.subset_go_validate_inputs <- function(
        go,
        ids
) {
    .assert_go_like(go = go)

    if (!is.character(ids) || length(ids) == 0L || anyNA(ids)) {
        rlang::abort(
            "`ids` must be a non-empty character vector of GO IDs.",
            arg = "ids"
        )
    }

    invisible(TRUE)
}


#' Extract GO IDs present in the graph
#'
#' @description
#' Retrieves the set of GO identifiers represented in a \code{GO} or
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
.extract_go_ids <- function(go) {
    ids <- go@terms$go_id
    ids <- ids[!is.na(ids)]
    unique(ids)
}


#' Check that seed GO IDs are present in the graph
#'
#' @description
#' Verifies that all supplied GO IDs occur in the current graph universe.
#' If any GO IDs are absent, an error is raised.
#'
#' @param ids \code{character()} vector of GO IDs to validate.
#' @param universe \code{character()} vector of GO IDs present in the
#'   current graph.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @noRd
.check_ids_in_graph <- function(
        ids,
        universe
) {
    bad <- setdiff(ids, universe)

    if (length(bad)) {
        msg <- sprintf(
            paste0(
                "Some `ids` are not present in this graph ",
                "(showing up to 10): %s"
            ),
            paste(utils::head(bad, 10L), collapse = ", ")
        )
        rlang::abort(message = msg, arg = "ids")
    }

    invisible(TRUE)
}


#' Expand GO seed terms to all descendants
#'
#' @description
#' Starting from one or more seed GO IDs, traverses the child adjacency
#' list to collect all descendants within the current graph universe.
#'
#' Seed terms are included in the returned set.
#'
#' @param seeds \code{character()} vector of seed GO IDs.
#' @param children Named list of child relationships from the GO graph.
#' @param universe \code{character()} vector of GO IDs present in the
#'   current graph.
#'
#' @return \code{character()} vector containing the seeds and all their
#'   descendants.
#'
#' @noRd
.expand_to_descendants <- function(
        seeds,
        children,
        universe
) {
    seen <- seeds
    queue <- seeds

    while (length(queue)) {
        node <- queue[[1L]]
        queue <- queue[-1L]

        kids <- children[[node]]
        if (is.null(kids) || !length(kids)) {
            next
        }

        kids <- kids[kids %in% universe]
        new <- setdiff(kids, seen)

        if (length(new)) {
            seen <- c(seen, new)
            queue <- c(queue, new)
        }
    }

    unique(seen)
}


#' Resolve retained GO IDs from seeds and mode
#'
#' @description
#' Determines the final set of retained GO IDs based on the descendant set
#' of the supplied seeds and the requested restriction mode.
#'
#' In \code{"keep"} mode, the descendant set is retained. In
#' \code{"exclude"} mode, the descendant set is removed from the graph
#' universe.
#'
#' @param universe \code{character()} vector of GO IDs present in the
#'   current graph.
#' @param seed_desc \code{character()} vector containing the seed GO IDs and
#'   all their descendants.
#' @param mode \code{character(1)} Restriction mode, either \code{"keep"} or
#'   \code{"exclude"}.
#'
#' @return \code{character()} vector of retained GO IDs.
#'
#' @noRd
.resolve_retained_ids <- function(
        universe,
        seed_desc,
        mode
) {
    if (identical(mode, "keep")) {
        return(universe[universe %in% seed_desc])
    }

    setdiff(universe, seed_desc)
}


#' Restrict GO term metadata to retained GO terms
#'
#' @description
#' Extracts the term metadata corresponding to the retained GO IDs after
#' graph restriction.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#' @param keep_ids \code{character()} vector of GO IDs to retain.
#'
#' @return A \code{data.frame} containing term metadata for
#'   \code{keep_ids}.
#'
#' @noRd
.restrict_terms <- function(
        go,
        keep_ids
) {
    terms_sub <- go@terms[match(keep_ids, go@terms$go_id), , drop = FALSE]
    terms_sub <- terms_sub[!is.na(terms_sub$go_id), , drop = FALSE]
    terms_sub
}
