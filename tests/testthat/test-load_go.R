# tests/testthat/test-load_go.R

.test_go_invariants <- function(obj, ont, ver) {
    testthat::expect_true(methods::is(obj, "GO"))
    testthat::expect_identical(obj@ontology, ont)
    testthat::expect_identical(obj@version, ver)

    testthat::expect_true(base::is.character(obj@ids))
    testthat::expect_true(base::is.data.frame(obj@terms))
    testthat::expect_true(base::nrow(obj@terms) > 0L)

    need <- c("go_id", "term", "ontology", "obsolete")
    testthat::expect_true(base::all(need %in% base::names(obj@terms)))
    testthat::expect_true(base::all(obj@terms$ontology == ont))
    testthat::expect_identical(obj@ids, obj@terms$go_id)

    testthat::expect_true(base::is.data.frame(obj@edges))
    testthat::expect_true(base::all(c("child", "parent") %in%
                                        base::names(obj@edges)))

    testthat::expect_true(base::is.list(obj@parents))
    testthat::expect_true(base::is.list(obj@children))
    testthat::expect_false(base::is.null(base::names(obj@parents)))
    testthat::expect_false(base::is.null(base::names(obj@children)))
    testthat::expect_identical(base::names(obj@parents), obj@ids)
    testthat::expect_identical(base::names(obj@children), obj@ids)

    testthat::expect_true(base::is.data.frame(obj@map))
    testthat::expect_true(base::all(c("go_id", "gene_id") %in%
                                        base::names(obj@map)))
    testthat::expect_identical(base::nrow(obj@map), 0L)

    if (base::nrow(obj@edges) > 0L) {
        keep <- obj@terms$go_id
        testthat::expect_true(base::all(obj@edges$child %in% keep))
        testthat::expect_true(base::all(obj@edges$parent %in% keep))
    }
}

testthat::test_that("GOcontext::load_go loads BP coherently by default", {
    testthat::skip_if_not_installed("GO.db")

    ver <- base::as.character(utils::packageVersion("GO.db"))
    obj <- GOcontext::load_go(ont = "BP")

    .test_go_invariants(obj, "BP", ver)
    testthat::expect_true(base::all(!obj@terms$obsolete))
})

testthat::test_that(
    "GOcontext::load_go include_obsolete=TRUE is accepted",
    {
        testthat::skip_if_not_installed("GO.db")

        ver <- base::as.character(utils::packageVersion("GO.db"))
        obj <- GOcontext::load_go(
            ont = "CC",
            include_obsolete = TRUE
        )

        .test_go_invariants(obj, "CC", ver)

        testthat::expect_true(
            base::all(c(TRUE, FALSE) %in% unique(obj@terms$obsolete)) ||
                base::all(!obj@terms$obsolete) ||
                base::all(obj@terms$obsolete)
        )
    }
)

testthat::test_that(
    "GOcontext::load_go rejects invalid include_obsolete values",
    {
        testthat::expect_error(
            GOcontext::load_go(
                ont = "CC",
                include_obsolete = NA
            ),
            regexp = "`include_obsolete` must be a logical\\(1\\)"
        )

        testthat::expect_error(
            GOcontext::load_go(
                ont = "CC",
                include_obsolete = "TRUE"
            ),
            regexp = "`include_obsolete` must be a logical\\(1\\)"
        )

        testthat::expect_error(
            GOcontext::load_go(
                ont = "CC",
                include_obsolete = c(TRUE, FALSE)
            ),
            regexp = "`include_obsolete` must be a logical\\(1\\)"
        )
    }
)

testthat::test_that(
    "GOcontext::load_go loads each ontology as a GO object",
    {
        testthat::skip_if_not_installed("GO.db")

        ver <- base::as.character(utils::packageVersion("GO.db"))

        bp <- GOcontext::load_go("BP")
        mf <- GOcontext::load_go("MF")
        cc <- GOcontext::load_go("CC")

        .test_go_invariants(bp, "BP", ver)
        .test_go_invariants(mf, "MF", ver)
        .test_go_invariants(cc, "CC", ver)
    }
)
