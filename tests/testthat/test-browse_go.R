.local_menu_choices <- function(choices) {
    i <- 0L

    function(...) {
        i <<- i + 1L
        if (i > base::length(choices)) {
            base::stop("Ran out of mocked menu choices.")
        }
        choices[[i]]
    }
}

.make_toy_go <- function() {
    ids <- c(
        "GO:0000001",
        "GO:0000002",
        "GO:0000003",
        "GO:0000004"
    )

    terms <- base::data.frame(
        go_id = ids,
        term = c("root", "mid", "leaf", "sib"),
        ontology = "CC",
        obsolete = c(FALSE, FALSE, FALSE, FALSE),
        stringsAsFactors = FALSE
    )

    edges <- base::data.frame(
        child = c("GO:0000002", "GO:0000003", "GO:0000004"),
        parent = c("GO:0000001", "GO:0000002", "GO:0000001"),
        stringsAsFactors = FALSE
    )

    parents <- list(
        "GO:0000001" = character(0),
        "GO:0000002" = "GO:0000001",
        "GO:0000003" = "GO:0000002",
        "GO:0000004" = "GO:0000001"
    )

    children <- list(
        "GO:0000001" = c("GO:0000002", "GO:0000004"),
        "GO:0000002" = "GO:0000003",
        "GO:0000003" = character(0),
        "GO:0000004" = character(0)
    )

    methods::new(
        "GO",
        ontology = "CC",
        version = "test-version",
        ids = ids,
        terms = terms,
        edges = edges,
        parents = parents,
        children = children,
        map = base::data.frame(
            go_id = character(0),
            gene_id = character(0),
            stringsAsFactors = FALSE
        )
    )
}

testthat::test_that(
    "GOcontext::browse_go rejects invalid input",
    {
        testthat::expect_error(
            GOcontext::browse_go("not_a_go_object"),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)

testthat::test_that(
    "GOcontext::browse_go returns NULL when menu is cancelled",
    {
        testthat::local_mocked_bindings(
            menu = .local_menu_choices(list(0L)),
            .package = "utils"
        )

        out <- GOcontext::browse_go(.go_cc)
        testthat::expect_null(out)
    }
)

testthat::test_that(
    "GOcontext::browse_go can finish immediately with no selection",
    {
        toy_go <- .make_toy_go()

        testthat::local_mocked_bindings(
            menu = .local_menu_choices(list(3L)),
            .package = "utils"
        )

        out <- GOcontext::browse_go(toy_go)

        testthat::expect_true(base::is.list(out))
        testthat::expect_identical(out$ids, character(0))
        testthat::expect_true(base::is.data.frame(out$summary))
        testthat::expect_identical(base::nrow(out$summary), 0L)
        testthat::expect_identical(
            base::names(out$summary),
            c("id", "term", "path")
        )
    }
)

testthat::test_that(
    "GOcontext::browse_go auto-adds a leaf term",
    {
        toy_go <- .make_toy_go()

        testthat::local_mocked_bindings(
            menu = .local_menu_choices(list(
                4L,  # select GO:0000002 from root view
                4L,  # select GO:0000003 leaf from mid view
                3L   # finish
            )),
            .package = "utils"
        )

        out <- GOcontext::browse_go(toy_go)

        testthat::expect_identical(out$ids, "GO:0000003")
        testthat::expect_identical(base::nrow(out$summary), 1L)
        testthat::expect_identical(out$summary$id, "GO:0000003")
        testthat::expect_identical(out$summary$term, "leaf")
    }
)

testthat::test_that(
    "GOcontext::browse_go adding an ancestor removes descendants",
    {
        toy_go <- .make_toy_go()

        testthat::local_mocked_bindings(
            menu = .local_menu_choices(list(
                4L,  # select GO:0000002 from root view
                4L,  # select GO:0000003 leaf from mid view, auto-add leaf
                4L,  # re-enter GO:0000002 from root view
                1L,  # add current node GO:0000002
                2L,  # go up to root view
                3L   # finish
            )),
            .package = "utils"
        )

        out <- GOcontext::browse_go(toy_go)

        testthat::expect_identical(out$ids, "GO:0000002")
        testthat::expect_identical(base::nrow(out$summary), 1L)
        testthat::expect_identical(out$summary$id, "GO:0000002")
        testthat::expect_identical(out$summary$term, "mid")
    }
)

testthat::test_that(
    "GOcontext::browse_go can add a non-root current node directly",
    {
        toy_go <- .make_toy_go()

        testthat::local_mocked_bindings(
            menu = .local_menu_choices(list(
                4L,  # select GO:0000002 from root view
                1L,  # add current node GO:0000002
                2L,  # go up
                3L   # finish
            )),
            .package = "utils"
        )

        out <- GOcontext::browse_go(toy_go)

        testthat::expect_identical(out$ids, "GO:0000002")
        testthat::expect_identical(base::nrow(out$summary), 1L)
        testthat::expect_identical(out$summary$term, "mid")
        testthat::expect_match(out$summary$path, "GO:0000001")
    }
)
