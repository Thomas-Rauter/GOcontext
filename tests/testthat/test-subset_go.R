testthat::test_that(
    "GOcontext::subset_go keep returns a coherent GOSubgraph",
    {
        seed <- .go_cc@terms$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = .go_cc,
            ids = seed,
            mode = "keep"
        )

        testthat::expect_true(methods::is(go_sub, "GOSubgraph"))
        testthat::expect_identical(go_sub@ontology, .go_cc@ontology)
        testthat::expect_identical(go_sub@version, .go_cc@version)
        testthat::expect_identical(go_sub@seed_ids, seed)
        testthat::expect_identical(go_sub@mode, "keep")

        testthat::expect_true(base::length(go_sub@keep_ids) >= 1L)
        testthat::expect_true(seed %in% go_sub@keep_ids)
        testthat::expect_identical(go_sub@keep_ids, go_sub@terms$go_id)

        testthat::expect_identical(
            base::sort(base::c(go_sub@keep_ids, go_sub@drop_ids)),
            base::sort(.go_cc@terms$go_id)
        )
        testthat::expect_identical(
            base::intersect(go_sub@keep_ids, go_sub@drop_ids),
            character(0)
        )

        testthat::expect_identical(
            base::names(go_sub@parents),
            go_sub@keep_ids
        )
        testthat::expect_identical(
            base::names(go_sub@children),
            go_sub@keep_ids
        )

        if (base::nrow(go_sub@edges) > 0L) {
            testthat::expect_true(
                base::all(go_sub@edges$child %in% go_sub@keep_ids)
            )
            testthat::expect_true(
                base::all(go_sub@edges$parent %in% go_sub@keep_ids)
            )
        }

        testthat::expect_identical(go_sub@map, .go_cc@map)
    }
)

testthat::test_that(
    "GOcontext::subset_go exclude removes the selected branch seed",
    {
        seed <- .go_cc@terms$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = .go_cc,
            ids = seed,
            mode = "exclude"
        )

        testthat::expect_true(methods::is(go_sub, "GOSubgraph"))
        testthat::expect_identical(go_sub@seed_ids, seed)
        testthat::expect_identical(go_sub@mode, "exclude")

        testthat::expect_false(seed %in% go_sub@keep_ids)
        testthat::expect_true(seed %in% go_sub@drop_ids)

        testthat::expect_identical(
            base::sort(base::c(go_sub@keep_ids, go_sub@drop_ids)),
            base::sort(.go_cc@terms$go_id)
        )
        testthat::expect_identical(
            base::intersect(go_sub@keep_ids, go_sub@drop_ids),
            character(0)
        )
    }
)

testthat::test_that(
    "GOcontext::subset_go removes duplicate seed IDs in output metadata",
    {
        seed <- .go_cc@terms$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = .go_cc,
            ids = base::c(seed, seed),
            mode = "keep"
        )

        testthat::expect_identical(go_sub@seed_ids, seed)
        testthat::expect_true(seed %in% go_sub@keep_ids)
    }
)

testthat::test_that(
    "GOcontext::subset_go works on an existing GOSubgraph",
    {
        seed1 <- .go_cc@terms$go_id[[1]]
        first_sub <- GOcontext::subset_go(
            go = .go_cc,
            ids = seed1,
            mode = "keep"
        )

        seed2 <- first_sub@terms$go_id[[1]]
        second_sub <- GOcontext::subset_go(
            go = first_sub,
            ids = seed2,
            mode = "keep"
        )

        testthat::expect_true(methods::is(second_sub, "GOSubgraph"))
        testthat::expect_true(
            base::all(second_sub@keep_ids %in% first_sub@terms$go_id)
        )
        testthat::expect_identical(
            base::names(second_sub@parents),
            second_sub@keep_ids
        )
        testthat::expect_identical(
            base::names(second_sub@children),
            second_sub@keep_ids
        )
    }
)

testthat::test_that(
    "GOcontext::subset_go preserves attached mapping when present",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        seed <- .go_cc_ecoli@terms$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = .go_cc_ecoli,
            ids = seed,
            mode = "keep"
        )

        testthat::expect_true(methods::is(go_sub, "GOSubgraph"))
        testthat::expect_identical(go_sub@map, .go_cc_ecoli@map)
    }
)

testthat::test_that(
    "GOcontext::subset_go rejects invalid GO objects",
    {
        testthat::expect_error(
            GOcontext::subset_go(
                go = "not_a_go_object",
                ids = "GO:0000001",
                mode = "keep"
            ),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)

testthat::test_that(
    "GOcontext::subset_go rejects invalid ids input",
    {
        testthat::expect_error(
            GOcontext::subset_go(
                go = .go_cc,
                ids = character(0),
                mode = "keep"
            ),
            regexp = "`ids` must be a non-empty character vector of GO IDs"
        )

        testthat::expect_error(
            GOcontext::subset_go(
                go = .go_cc,
                ids = NA_character_,
                mode = "keep"
            ),
            regexp = "`ids` must be a non-empty character vector of GO IDs"
        )

        testthat::expect_error(
            GOcontext::subset_go(
                go = .go_cc,
                ids = 1,
                mode = "keep"
            ),
            regexp = "`ids` must be a non-empty character vector of GO IDs"
        )
    }
)

testthat::test_that(
    "GOcontext::subset_go rejects IDs not present in the graph",
    {
        testthat::expect_error(
            GOcontext::subset_go(
                go = .go_cc,
                ids = "GO:9999999",
                mode = "keep"
            ),
            regexp = "Some `ids` are not present in this graph"
        )
    }
)
