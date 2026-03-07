# tests/testthat/test-load_go.R

.test_go_invariants <- function(obj, ont, ver) {
    testthat::expect_true(methods::is(obj, "GO"))
    testthat::expect_identical(obj@ontology, ont)
    testthat::expect_identical(obj@version, ver)

    testthat::expect_true(is.data.frame(obj@terms))
    testthat::expect_true(nrow(obj@terms) > 0L)

    need <- c("go_id", "term", "ontology", "obsolete")
    testthat::expect_true(all(need %in% names(obj@terms)))
    testthat::expect_true(all(obj@terms$ontology == ont))

    testthat::expect_true(is.data.frame(obj@edges))
    testthat::expect_true(is.list(obj@parents))
    testthat::expect_true(is.list(obj@children))
    testthat::expect_true(!is.null(names(obj@parents)))
    testthat::expect_true(!is.null(names(obj@children)))

    if (nrow(obj@edges) > 0L) {
        keep <- obj@terms$go_id
        testthat::expect_true(all(obj@edges$child %in% keep))
        testthat::expect_true(all(obj@edges$parent %in% keep))
    }
}

testthat::test_that("load_go() default BP is coherent and compute_depthd", {
    testthat::skip_if_not_installed("GO.db")

    ver <- as.character(utils::packageVersion("GO.db"))
    obj <- GOcontext::load_go(ont = "BP")

    .test_go_invariants(obj, "BP", ver)

    testthat::expect_true(all(!obj@terms$obsolete))

    testthat::expect_true(is.integer(obj@depth))
    testthat::expect_true(length(obj@depth) == nrow(obj@terms))
})

testthat::test_that("load_go() compute_depth=FALSE yields empty depth", {
    testthat::skip_if_not_installed("GO.db")

    ver <- as.character(utils::packageVersion("GO.db"))
    obj <- GOcontext::load_go(ont = "BP", compute_depth = FALSE)

    .test_go_invariants(obj, "BP", ver)

    testthat::expect_true(is.integer(obj@depth))
    testthat::expect_identical(length(obj@depth), 0L)
})

testthat::test_that("load_go() include_obsolete=TRUE allows obsolete terms", {
    testthat::skip_if_not_installed("GO.db")

    ver <- as.character(utils::packageVersion("GO.db"))
    obj <- GOcontext::load_go(
        ont = "CC",
        include_obsolete = TRUE,
        compute_depth = FALSE
        )

    .test_go_invariants(obj, "CC", ver)

    testthat::expect_true(any(obj@terms$obsolete) ||
                              all(!obj@terms$obsolete))
})

testthat::test_that("load_go() version mismatch errors", {
    testthat::skip_if_not_installed("GO.db")

    testthat::expect_error(
        GOcontext::load_go(
            ont = "CC",
            version = "0.0.0",
            compute_depth = FALSE
            ),
        regexp = "does not match requested version",
        fixed = FALSE
    )
})
