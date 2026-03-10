.make_empty_go <- function() {
    terms <- base::data.frame(
        go_id = character(0),
        term = character(0),
        ontology = character(0),
        obsolete = logical(0),
        stringsAsFactors = FALSE
    )

    methods::new(
        "GO",
        ontology = "CC",
        version = "test-version",
        ids = character(0),
        terms = terms,
        edges = base::data.frame(
            child = character(0),
            parent = character(0),
            stringsAsFactors = FALSE
        ),
        parents = stats::setNames(list(), character(0)),
        children = stats::setNames(list(), character(0)),
        map = base::data.frame(
            go_id = character(0),
            gene_id = character(0),
            stringsAsFactors = FALSE
        )
    )
}

testthat::test_that(
    "GOcontext:::.empty_go_map returns the expected empty mapping table",
    {
        out <- GOcontext:::.empty_go_map()

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("go_id", "gene_id")
        )
        testthat::expect_identical(base::nrow(out), 0L)
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_validate_inputs accepts valid inputs",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        testthat::expect_invisible(
            GOcontext:::.attach_org_validate_inputs(
                go = .go_cc,
                OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                keytype = "ENTREZID"
            )
        )
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_validate_inputs rejects invalid OrgDb",
    {
        testthat::expect_error(
            GOcontext:::.attach_org_validate_inputs(
                go = .go_cc,
                OrgDb = "not_an_orgdb",
                keytype = "ENTREZID"
            ),
            regexp = "`OrgDb` must be an OrgDb object"
        )
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_validate_inputs rejects invalid keytype",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        org_db <- org.EcK12.eg.db::org.EcK12.eg.db

        testthat::expect_error(
            GOcontext:::.attach_org_validate_inputs(
                go = .go_cc,
                OrgDb = org_db,
                keytype = NA_character_
            ),
            regexp = "`keytype` must be a non-NA character\\(1\\)"
        )

        testthat::expect_error(
            GOcontext:::.attach_org_validate_inputs(
                go = .go_cc,
                OrgDb = org_db,
                keytype = ""
            ),
            regexp = "`keytype` must not be empty"
        )

        testthat::expect_error(
            GOcontext:::.attach_org_validate_inputs(
                go = .go_cc,
                OrgDb = org_db,
                keytype = "NOT_A_KEYTYPE"
            ),
            regexp = "is not available in `OrgDb`"
        )
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_build_map returns a clean mapping table",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        out <- GOcontext:::.attach_org_build_map(
            go = .go_cc,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("go_id", "gene_id")
        )
        testthat::expect_true(base::nrow(out) > 0L)
        testthat::expect_false(base::anyNA(out$go_id))
        testthat::expect_false(base::anyNA(out$gene_id))
        testthat::expect_identical(base::anyDuplicated(out), 0L)
        testthat::expect_true(
            base::all(out$go_id %in% .go_cc@terms$go_id)
        )
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_build_map returns empty map for empty GO",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_empty <- .make_empty_go()

        out <- GOcontext:::.attach_org_build_map(
            go = go_empty,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        testthat::expect_true(base::is.data.frame(out))
        testthat::expect_identical(
            base::names(out),
            c("go_id", "gene_id")
        )
        testthat::expect_identical(base::nrow(out), 0L)
    }
)

testthat::test_that(
    "GOcontext:::.attach_org_set_map writes the map into the object",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        map <- GOcontext:::.attach_org_build_map(
            go = .go_cc,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        go_out <- GOcontext:::.attach_org_set_map(
            go = .go_cc,
            map = map
        )

        testthat::expect_true(methods::is(go_out, "GO"))
        testthat::expect_identical(go_out@map, map)
        testthat::expect_identical(go_out@ontology, .go_cc@ontology)
        testthat::expect_identical(go_out@version, .go_cc@version)
        testthat::expect_identical(go_out@ids, .go_cc@ids)
        testthat::expect_identical(go_out@terms, .go_cc@terms)
    }
)
