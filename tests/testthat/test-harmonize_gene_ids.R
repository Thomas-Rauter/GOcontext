testthat::test_that(
    "GOcontext::harmonize_gene_ids harmonizes valid identifiers",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        genes <- c("b0002", "thrA", "b0003")

        res <- GOcontext::harmonize_gene_ids(
            genes = genes,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            to_keytype = "ENTREZID",
            from_keytypes = c("SYMBOL", "ALIAS")
        )

        testthat::expect_true(base::is.data.frame(res))
        testthat::expect_identical(base::nrow(res), length(genes))
        testthat::expect_true(
            base::all(
                c(
                    "input",
                    "matches_target",
                    "from_match_status",
                    "to_match_status",
                    "harmonized"
                ) %in% base::colnames(res)
            )
        )
        testthat::expect_identical(res$input, genes)
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids preserves duplicates and NA rows",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        genes <- c("b0002", "b0002", NA_character_, "not_a_gene")

        res <- GOcontext::harmonize_gene_ids(
            genes = genes,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            to_keytype = "ENTREZID",
            from_keytypes = c("SYMBOL", "ALIAS")
        )

        testthat::expect_identical(base::nrow(res), length(genes))
        testthat::expect_identical(res$input, genes)

        testthat::expect_true(base::is.na(res$harmonized[[3]]))
        testthat::expect_true(base::is.na(res$from_match_status[[3]]))
        testthat::expect_true(base::is.na(res$to_match_status[[3]]))
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids identifies direct target matches",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        genes <- c("944742")

        res <- GOcontext::harmonize_gene_ids(
            genes = genes,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            to_keytype = "ENTREZID",
            from_keytypes = c("SYMBOL", "ALIAS")
        )

        testthat::expect_true(res$matches_target[[1]])
        testthat::expect_true(base::is.na(res$from_match_status[[1]]))
        testthat::expect_true(base::is.na(res$to_match_status[[1]]))
        testthat::expect_identical(res$harmonized[[1]], genes[[1]])
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids reports unresolved identifiers",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        genes <- c("not_a_gene")

        res <- GOcontext::harmonize_gene_ids(
            genes = genes,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            to_keytype = "ENTREZID",
            from_keytypes = c("SYMBOL", "ALIAS")
        )

        testthat::expect_identical(res$from_match_status[[1]], "none")
        testthat::expect_true(base::is.na(res$harmonized[[1]]))
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids rejects invalid OrgDb objects",
    {
        testthat::expect_error(
            GOcontext::harmonize_gene_ids(
                genes = c("b0002"),
                OrgDb = "not_an_orgdb",
                to_keytype = "ENTREZID",
                from_keytypes = c("SYMBOL")
            ),
            regexp = "`OrgDb` must be an OrgDb object"
        )
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids rejects invalid genes input",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::harmonize_gene_ids(
                genes = 123,
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                to_keytype = "ENTREZID",
                from_keytypes = c("SYMBOL")
            ),
            regexp = "`genes` must be a character vector"
        )
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids rejects invalid to_keytype values",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::harmonize_gene_ids(
                genes = c("b0002"),
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                to_keytype = NA_character_,
                from_keytypes = c("SYMBOL")
            ),
            regexp = "`to_keytype` must be a non-empty character"
        )

        testthat::expect_error(
            GOcontext::harmonize_gene_ids(
                genes = c("b0002"),
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                to_keytype = "NOT_A_KEYTYPE",
                from_keytypes = c("SYMBOL")
            ),
            regexp = "is not available in `OrgDb`"
        )
    }
)

testthat::test_that(
    "GOcontext::harmonize_gene_ids rejects invalid from_keytypes",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_error(
            GOcontext::harmonize_gene_ids(
                genes = c("b0002"),
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                to_keytype = "ENTREZID",
                from_keytypes = c("SYMBOL", NA_character_)
            ),
            regexp = "`from_keytypes` must be NULL or a character vector"
        )
    }
)
