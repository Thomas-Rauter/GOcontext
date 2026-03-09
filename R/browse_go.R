#' Interactively browse a GO graph
#'
#' @description
#' Provides an interactive console browser for exploring the Gene Ontology
#' hierarchy. The user can navigate through parent–child relationships,
#' add nodes of interest, and return the selected GO IDs.
#'
#' @param go A \code{GO} or \code{GOSubgraph} object.
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item \code{ids}: \code{character()} vector of selected GO IDs.
#'   \item \code{summary}: \code{data.frame} summarizing the selected terms.
#' }
#'
#' @export
browse_go <- function(go) {
    .assert_go_like(go = go)

    children <- go@children
    parents  <- go@parents
    terms    <- go@terms

    roots <- names(parents)[lengths(parents) == 0]

    # start at children of ontology roots
    current <- unique(unlist(children[roots], use.names = FALSE))
    current_node <- roots[1]

    stack <- list()
    selected <- character()

    repeat {
        labels <- terms$term[match(current, terms$go_id)]
        menu_entries <- paste0(current, "  ", labels)

        choices <- c(
            "[add current node]",
            "[go up]",
            "[finish]",
            menu_entries
        )

        title <- paste("GO browser:", current_node)

        sel <- utils::menu(choices, title = title)
        if (sel == 0) return(invisible(NULL))
        if (sel == 3) break
        if (sel == 2) {
            if (length(stack) > 0) {
                state <- stack[[length(stack)]]
                stack <- stack[-length(stack)]
                current <- state$current
                current_node <- state$current_node
            }

            next
        }

        if (sel == 1) {

            if (!is.na(current_node) && !(current_node %in% roots)) {

                covered_by <- .is_descendant_of_any(
                    current_node,
                    selected,
                    parents
                    )
                if (!is.na(covered_by)) {
                    message(
                        "Not added. Already covered by ancestor: ",
                        covered_by
                        )
                    next
                }

                covered_desc <- .is_ancestor_of_any(
                    current_node,
                    selected,
                    parents
                    )
                if (length(covered_desc)) {
                    selected <- setdiff(selected, covered_desc)
                    message(
                        "Removed covered selections: ",
                        paste(covered_desc, collapse = ", ")
                        )
                }

                selected <- unique(c(selected, current_node))
                message("Added: ", current_node)
            }

            next
        }

        idx <- sel - 3
        node <- current[idx]

        stack[[length(stack) + 1]] <- list(
            current = current,
            current_node = current_node
        )

        current_node <- node
        current <- children[[node]]

        if (!length(current)) {
            message(
                "Term ",
                node,
                " has no children in this graph. ",
                "Added to the list of selected terms."
                )
            selected <- unique(c(selected, node))

            state <- stack[[length(stack)]]
            stack <- stack[-length(stack)]

            current <- state$current
            current_node <- state$current_node
        }
    }

    summary <- .create_summary(
        selected = selected,
        terms = terms,
        parents = parents
    )

    list(
        ids = selected,
        summary = summary
    )
}


# Level 1 helpers --------------------------------------------------------------


#' Check whether a GO term is an ancestor of another
#'
#' @description
#' Determines whether GO term \code{a} is an ancestor of GO term \code{b}
#' by traversing parent relationships in the GO graph.
#'
#' Because the Gene Ontology is a directed acyclic graph (DAG), terms may
#' have multiple parents. This function performs a breadth-first traversal
#' of parent relationships starting from \code{b} and returns \code{TRUE}
#' if \code{a} occurs on any ancestor path.
#'
#' The comparison is strict: if \code{a == b}, the function returns
#' \code{FALSE}.
#'
#' @param a \code{character(1)} Candidate ancestor GO ID.
#' @param b \code{character(1)} Candidate descendant GO ID.
#' @param parents Named list of parent relationships from the GO graph.
#'
#' @return \code{logical(1)} indicating whether \code{a} is an ancestor of
#'   \code{b}.
#'
#' @noRd
.is_ancestor_of <- function(
        a,
        b,
        parents
) {
    if (identical(a, b)) {
        return(FALSE)
    }

    seen <- character(0)
    queue <- b

    while (length(queue)) {
        cur <- queue[[1L]]
        queue <- queue[-1L]

        pr <- parents[[cur]]
        if (!length(pr)) {
            next
        }

        if (a %in% pr) {
            return(TRUE)
        }

        new <- setdiff(pr, seen)
        if (length(new)) {
            seen <- c(seen, new)
            queue <- c(queue, new)
        }
    }

    FALSE
}


#' Check whether a GO term is a descendant of any selected term
#'
#' @description
#' Tests whether a GO term is a descendant of any GO ID in a set of
#' selected terms. Used to avoid adding nodes that are already covered
#' by a previously selected ancestor.
#'
#' @param x \code{character(1)} GO ID to test.
#' @param sel \code{character()} Vector of selected GO IDs.
#' @param parents Named list of parent relationships from the GO graph.
#'
#' @return The GO ID of the covering ancestor if one exists, otherwise
#'   \code{NA_character_}.
#'
#' @noRd
.is_descendant_of_any <- function(
        x,
        sel,
        parents
) {
    for (s in sel) {
        if (.is_ancestor_of(s, x, parents)) return(s)
    }
    NA_character_
}


#' Find selected GO terms that are descendants of a given term
#'
#' @description
#' Identifies selected GO IDs that are descendants of the supplied GO
#' term. Used when adding a new node to remove selections that become
#' redundant because they are covered by the new ancestor.
#'
#' @param x \code{character(1)} Candidate ancestor GO ID.
#' @param sel \code{character()} Vector of selected GO IDs.
#' @param parents Named list of parent relationships from the GO graph.
#'
#' @return \code{character()} vector of selected GO IDs that are
#'   descendants of \code{x}.
#'
#' @noRd
.is_ancestor_of_any <- function(x, sel, parents) {
    sel[vapply(sel, function(s) .is_ancestor_of(x, s, parents), logical(1))]
}


#' Build a summary table for selected GO terms
#'
#' @description
#' Constructs a summary table describing the GO terms selected during
#' interactive browsing. The table contains GO IDs, term labels, and a
#' hierarchical path representation for each term.
#'
#' @param selected \code{character()} Vector of selected GO IDs.
#' @param terms GO term metadata table containing at least columns
#'   \code{go_id} and \code{term}.
#' @param parents Named list of parent relationships from the GO graph.
#'
#' @return A \code{data.frame} with columns \code{id}, \code{term}, and
#'   \code{path}.
#'
#' @noRd
.create_summary <- function(
        selected,
        terms,
        parents
) {
    data.frame(
        id = selected,
        term = terms$term[match(selected, terms$go_id)],
        path = vapply(
            selected,
            FUN = .go_path_string,
            FUN.VALUE = character(1),
            parents = parents
        ),
        stringsAsFactors = FALSE
    )
}


# Level 2 helpers --------------------------------------------------------------


#' Construct a hierarchical path string for a GO term
#'
#' @description
#' Builds a textual representation of one hierarchical path from a GO
#' term to an ontology root by traversing parent relationships.
#'
#' Because the Gene Ontology is a directed acyclic graph (DAG), terms
#' may have multiple parents. This function follows one parent path to
#' construct a representative rootward path string.
#'
#' The path is returned as a string with GO IDs separated by
#' \code{" > "}.
#'
#' @param id \code{character(1)} GO ID for which to compute the path.
#' @param parents Named list of parent relationships from the GO graph.
#'
#' @return \code{character(1)} containing the hierarchical path.
#'
#' @noRd
.go_path_string <- function(
        id,
        parents
) {
    p <- id
    cur <- id

    repeat {
        pr <- parents[[cur]]
        if (length(pr) == 0L) {
            break
        }
        cur <- pr[1]
        p <- c(cur, p)
    }

    paste(p, collapse = " > ")
}
