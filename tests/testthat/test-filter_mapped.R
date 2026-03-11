testthat::test_that(
    "GOcontext::filter_mapped returns a mapped GOSubgraph",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_sub <- GOcontext::filter_mapped(go_cc_ecoli)

        expected_seed <- base::intersect(
            base::unique(go_cc_ecoli@map$go_id),
            go_cc_ecoli@terms$go_id
        )

        testthat::expect_true(methods::is(go_sub, "GOSubgraph"))
        testthat::expect_identical(go_sub@ontology, go_cc_ecoli@ontology)
        testthat::expect_identical(go_sub@version, go_cc_ecoli@version)
        testthat::expect_identical(
            base::sort(go_sub@seed_ids),
            base::sort(expected_seed)
        )
        testthat::expect_identical(go_sub@mode, "keep")

        testthat::expect_true(base::is.character(go_sub@keep_ids))
        testthat::expect_true(base::is.character(go_sub@drop_ids))
        testthat::expect_true(base::is.data.frame(go_sub@terms))
        testthat::expect_true(base::is.data.frame(go_sub@edges))
        testthat::expect_true(base::is.list(go_sub@parents))
        testthat::expect_true(base::is.list(go_sub@children))

        testthat::expect_identical(go_sub@keep_ids, go_sub@terms$go_id)
        testthat::expect_identical(
            base::sort(base::c(go_sub@keep_ids, go_sub@drop_ids)),
            base::sort(go_cc_ecoli@terms$go_id)
        )
        testthat::expect_identical(
            base::intersect(go_sub@keep_ids, go_sub@drop_ids),
            character(0)
        )

        testthat::expect_true(
            base::all(go_sub@seed_ids %in% go_sub@keep_ids)
        )
        testthat::expect_true(
            base::all(go_sub@seed_ids %in% go_cc_ecoli@map$go_id)
        )

        testthat::expect_identical(go_sub@map, go_cc_ecoli@map)

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
    }
)

testthat::test_that(
    "GOcontext::filter_mapped keeps mapped terms and their ancestors",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_sub <- GOcontext::filter_mapped(go_cc_ecoli)

        expected_seed <- base::intersect(
            base::unique(go_cc_ecoli@map$go_id),
            go_cc_ecoli@terms$go_id
        )

        testthat::expect_identical(
            base::sort(go_sub@seed_ids),
            base::sort(expected_seed)
        )

        testthat::expect_true(
            base::all(expected_seed %in% go_sub@keep_ids)
        )
        testthat::expect_true(
            length(go_sub@keep_ids) >= length(expected_seed)
        )
    }
)

testthat::test_that(
    "GOcontext::filter_mapped works on an already restricted graph",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        mapped_id <- go_cc_ecoli@map$go_id[[1]]

        go_small <- GOcontext::subset_go(
            go = go_cc_ecoli,
            ids = mapped_id,
            mode = "keep"
        )

        go_small_mapped <- GOcontext::filter_mapped(go_small)

        expected_seed <- base::intersect(
            base::unique(go_small@map$go_id),
            go_small@terms$go_id
        )

        testthat::expect_true(methods::is(go_small_mapped, "GOSubgraph"))
        testthat::expect_identical(
            go_small_mapped@ontology,
            go_small@ontology
        )
        testthat::expect_identical(
            go_small_mapped@version,
            go_small@version
        )
        testthat::expect_identical(
            base::sort(go_small_mapped@seed_ids),
            base::sort(expected_seed)
        )
        testthat::expect_identical(go_small_mapped@mode, "keep")

        testthat::expect_true(
            base::all(go_small_mapped@keep_ids %in% go_small@terms$go_id)
        )
        testthat::expect_true(
            base::all(go_small_mapped@seed_ids %in% go_small@map$go_id)
        )
        testthat::expect_true(
            base::all(go_small_mapped@seed_ids %in% go_small_mapped@keep_ids)
        )
        testthat::expect_identical(go_small_mapped@map, go_small@map)
    }
)

testthat::test_that(
    "GOcontext::filter_mapped rejects GO objects without mapping",
    {
        testthat::expect_error(
            GOcontext::filter_mapped(go_cc),
            regexp = "must have a non-empty attached organism mapping"
        )
    }
)

testthat::test_that(
    "GOcontext::filter_mapped rejects invalid input objects",
    {
        testthat::expect_error(
            GOcontext::filter_mapped("not_a_go_object"),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)
