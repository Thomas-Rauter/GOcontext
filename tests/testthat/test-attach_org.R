testthat::test_that(
    "GOcontext::attach_org attaches an E. coli GO-to-gene map",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_out <- GOcontext::attach_org(
            go = go_cc,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        testthat::expect_true(methods::is(go_out, "GO"))
        testthat::expect_true(base::is.data.frame(go_out@map))
        testthat::expect_true(
            base::all(c("go_id", "gene_id") %in% base::colnames(go_out@map))
        )
        testthat::expect_gt(base::nrow(go_out@map), 0L)
        testthat::expect_false(base::anyNA(go_out@map$go_id))
        testthat::expect_false(base::anyNA(go_out@map$gene_id))
        testthat::expect_identical(base::anyDuplicated(go_out@map), 0L)
        testthat::expect_true(
            base::all(go_out@map$go_id %in% go_out@terms$go_id)
        )

        testthat::expect_identical(go_out@ontology, go_cc@ontology)
        testthat::expect_identical(go_out@version, go_cc@version)
        testthat::expect_identical(go_out@ids, go_cc@ids)
        testthat::expect_identical(go_out@terms, go_cc@terms)
        testthat::expect_identical(go_out@edges, go_cc@edges)
        testthat::expect_identical(go_out@parents, go_cc@parents)
        testthat::expect_identical(go_out@children, go_cc@children)
    }
)

testthat::test_that(
    "shared E. coli-attached fixture is available for downstream tests",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_true(base::exists("go_cc_ecoli", inherits = TRUE))
        testthat::expect_true(methods::is(go_cc_ecoli, "GO"))
        testthat::expect_true(base::is.data.frame(go_cc_ecoli@map))
        testthat::expect_gt(base::nrow(go_cc_ecoli@map), 0L)
        testthat::expect_true(
            base::all(go_cc_ecoli@map$go_id %in% go_cc_ecoli@terms$go_id)
        )
    }
)

testthat::test_that(
    "GOcontext::attach_org works on a GOSubgraph with E. coli OrgDb",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        mapped_id <- go_cc_ecoli@map$go_id[[1]]

        go_sub <- GOcontext::subset_go(
            go = go_cc,
            ids = mapped_id,
            mode = "keep"
        )

        go_sub_out <- GOcontext::attach_org(
            go = go_sub,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        testthat::expect_true(methods::is(go_sub_out, "GOSubgraph"))
        testthat::expect_true(
            base::all(go_sub_out@map$go_id %in% go_sub_out@terms$go_id)
        )
        testthat::expect_identical(go_sub_out@keep_ids, go_sub@keep_ids)
        testthat::expect_identical(go_sub_out@drop_ids, go_sub@drop_ids)
        testthat::expect_identical(go_sub_out@seed_ids, go_sub@seed_ids)
        testthat::expect_identical(go_sub_out@mode, go_sub@mode)
        testthat::expect_identical(go_sub_out@terms, go_sub@terms)
        testthat::expect_identical(go_sub_out@edges, go_sub@edges)
        testthat::expect_identical(go_sub_out@parents, go_sub@parents)
        testthat::expect_identical(go_sub_out@children, go_sub@children)
    }
)

testthat::test_that(
    "GOcontext::attach_org rejects an invalid GO object",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::attach_org(
                go = "not_a_go_object",
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                keytype = "ENTREZID"
            ),
            regexp = "must be a GO or GOSubgraph object"
        )
    }
)

testthat::test_that(
    "GOcontext::attach_org rejects an invalid OrgDb object",
    {
        testthat::expect_error(
            GOcontext::attach_org(
                go = go_cc,
                OrgDb = "not_an_orgdb",
                keytype = "ENTREZID"
            ),
            regexp = "`OrgDb` must be an OrgDb object"
        )
    }
)

testthat::test_that(
    "GOcontext::attach_org rejects invalid keytype values",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::attach_org(
                go = go_cc,
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                keytype = NA_character_
            ),
            regexp = "`keytype` must be a non-NA character\\(1\\)"
        )

        testthat::expect_error(
            GOcontext::attach_org(
                go = go_cc,
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                keytype = ""
            ),
            regexp = "`keytype` must not be empty"
        )

        testthat::expect_error(
            GOcontext::attach_org(
                go = go_cc,
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                keytype = "NOT_A_KEYTYPE"
            ),
            regexp = "is not available in `OrgDb`"
        )
    }
)

testthat::test_that(
    "GOcontext::attach_org trims surrounding whitespace in keytype",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_out <- GOcontext::attach_org(
            go = go_cc,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = " ENTREZID "
        )

        testthat::expect_true(base::is.data.frame(go_out@map))
        testthat::expect_gt(base::nrow(go_out@map), 0L)
        testthat::expect_true(
            base::all(go_out@map$go_id %in% go_out@terms$go_id)
        )
    }
)
