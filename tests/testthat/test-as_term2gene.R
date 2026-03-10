testthat::test_that(
    "GOcontext::as_term2gene returns a size-filtered TERM2GENE table for GO",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        out <- GOcontext::as_term2gene(go_cc_ecoli)

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "gene")
        )
        testthat::expect_true(base::nrow(out) > 0L)

        testthat::expect_true(base::is.character(out$term))
        testthat::expect_true(base::is.character(out$gene))
        testthat::expect_false(base::anyNA(out$term))
        testthat::expect_false(base::anyNA(out$gene))
        testthat::expect_identical(base::anyDuplicated(out), 0L)

        testthat::expect_true(
            base::all(out$term %in% go_cc_ecoli@terms$go_id)
        )
        testthat::expect_true(
            base::all(out$term %in% go_cc_ecoli@map$go_id)
        )
        testthat::expect_true(
            base::all(out$gene %in% go_cc_ecoli@map$gene_id)
        )

        term_sizes <- table(go_cc_ecoli@map$go_id)
        keep_terms <- names(term_sizes)[
            term_sizes >= 10L & term_sizes <= 500L
        ]

        testthat::expect_true(base::all(out$term %in% keep_terms))
    }
)

testthat::test_that(
    "GOcontext::as_term2gene returns a size-filtered TERM2GENE table for GOSubgraph",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        seed <- go_cc_ecoli@map$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = go_cc_ecoli,
            ids = seed,
            mode = "keep"
        )

        out <- GOcontext::as_term2gene(
            go_sub,
            minGSSize = 1L,
            maxGSSize = 500L
        )

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "gene")
        )
        testthat::expect_true(base::nrow(out) > 0L)
        testthat::expect_identical(base::anyDuplicated(out), 0L)

        testthat::expect_true(
            base::all(out$term %in% go_sub@terms$go_id)
        )
        testthat::expect_true(
            base::all(out$term %in% go_sub@map$go_id)
        )

        term_sizes <- table(
            go_sub@map$go_id[go_sub@map$go_id %in% go_sub@terms$go_id]
        )
        keep_terms <- names(term_sizes)[
            term_sizes >= 1L & term_sizes <= 500L
        ]

        testthat::expect_true(base::all(out$term %in% keep_terms))
    }
)

testthat::test_that(
    "GOcontext::as_term2gene respects custom gene set size thresholds",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        out <- GOcontext::as_term2gene(
            go_cc_ecoli,
            minGSSize = 1L,
            maxGSSize = 5L
        )

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "gene")
        )

        if (base::nrow(out) > 0L) {
            term_sizes <- table(out$term)
            testthat::expect_true(base::all(term_sizes >= 1L))
            testthat::expect_true(base::all(term_sizes <= 5L))
        }
    }
)

testthat::test_that(
    "GOcontext::as_term2gene rejects invalid size thresholds",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::as_term2gene(
                go_cc_ecoli,
                minGSSize = 0
            ),
            regexp = "`minGSSize` must be a positive integer"
        )

        testthat::expect_error(
            GOcontext::as_term2gene(
                go_cc_ecoli,
                minGSSize = 10L,
                maxGSSize = 5L
            ),
            regexp = "`maxGSSize` must be >= `minGSSize`"
        )
    }
)

testthat::test_that(
    "GOcontext::as_term2gene rejects GO objects without mapping",
    {
        testthat::expect_error(
            GOcontext::as_term2gene(go_cc),
            regexp = "must have a non-empty attached organism mapping"
        )
    }
)

testthat::test_that(
    "GOcontext::as_term2gene rejects invalid input objects",
    {
        testthat::expect_error(
            GOcontext::as_term2gene("not_a_go_object"),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)

testthat::test_that(
    "GOcontext::as_term2gene drops missing and duplicate mappings",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go <- go_cc_ecoli
        old_map <- go@map

        on.exit(
            {
                go@map <- old_map
            },
            add = TRUE
        )

        clean_out <- GOcontext::as_term2gene(go)

        extra_na_go <- old_map[1, , drop = FALSE]
        extra_na_go$go_id <- NA_character_

        extra_na_gene <- old_map[1, , drop = FALSE]
        extra_na_gene$gene_id <- NA_character_

        dup <- old_map[1, , drop = FALSE]

        go@map <- base::rbind(old_map, extra_na_go, extra_na_gene, dup)

        out <- GOcontext::as_term2gene(go)

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("term", "gene")
        )
        testthat::expect_false(base::anyNA(out$term))
        testthat::expect_false(base::anyNA(out$gene))
        testthat::expect_identical(base::anyDuplicated(out), 0L)
        testthat::expect_identical(out, clean_out)
    }
)

testthat::test_that(
    "GOcontext::as_term2gene restricts mappings to the current graph",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        seed <- go_cc_ecoli@map$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = go_cc_ecoli,
            ids = seed,
            mode = "keep"
        )

        out <- GOcontext::as_term2gene(go_sub)

        testthat::expect_true(
            base::all(out$term %in% go_sub@terms$go_id)
        )
        testthat::expect_false(
            base::any(out$term %in% go_sub@drop_ids)
        )
    }
)
