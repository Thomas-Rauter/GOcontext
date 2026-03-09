testthat::test_that(
    "GOcontext::as_term2name returns a TERM2NAME table for GO",
    {
        out <- GOcontext::as_term2name(.go_cc)

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "name")
        )
        testthat::expect_identical(
            base::nrow(out),
            base::nrow(.go_cc@terms)
        )

        testthat::expect_true(base::is.character(out$term))
        testthat::expect_true(base::is.character(out$name))
        testthat::expect_false(base::anyNA(out$term))
        testthat::expect_false(base::anyNA(out$name))
        testthat::expect_identical(base::anyDuplicated(out$term), 0L)

        testthat::expect_identical(
            out$term,
            .go_cc@terms$go_id
        )
        testthat::expect_identical(
            out$name,
            .go_cc@terms$term
        )
    }
)

testthat::test_that(
    "GOcontext::as_term2name returns a TERM2NAME table for GOSubgraph",
    {
        seed <- .go_cc@terms$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = .go_cc,
            ids = seed,
            mode = "keep"
        )

        out <- GOcontext::as_term2name(go_sub)

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "name")
        )
        testthat::expect_identical(
            base::nrow(out),
            base::nrow(go_sub@terms)
        )

        testthat::expect_identical(
            out$term,
            go_sub@terms$go_id
        )
        testthat::expect_identical(
            out$name,
            go_sub@terms$term
        )
        testthat::expect_identical(base::anyDuplicated(out$term), 0L)
    }
)

testthat::test_that(
    "GOcontext::as_term2name rejects invalid input objects",
    {
        testthat::expect_error(
            GOcontext::as_term2name("not_a_go_object"),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)

testthat::test_that(
    "GOcontext::as_term2name drops missing and duplicate term rows",
    {
        go <- .go_cc
        old_terms <- go@terms
        old_ids <- go@ids

        on.exit(
            {
                go@terms <- old_terms
                go@ids <- old_ids
            },
            add = TRUE
        )

        extra <- old_terms[1, , drop = FALSE]
        extra$term <- NA_character_

        dup <- old_terms[1, , drop = FALSE]

        go@terms <- base::rbind(old_terms, extra, dup)
        go@ids <- go@terms$go_id

        out <- GOcontext::as_term2name(go)

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "name")
        )
        testthat::expect_false(base::anyNA(out$term))
        testthat::expect_false(base::anyNA(out$name))
        testthat::expect_identical(base::anyDuplicated(out$term), 0L)
        testthat::expect_identical(
            base::nrow(out),
            base::nrow(old_terms)
        )
    }
)
