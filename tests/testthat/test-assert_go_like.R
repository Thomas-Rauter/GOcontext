# tests/testthat/test-assert_go_like.R

testthat::test_that(".assert_go_like() errors on invalid ontology slot", {
    go <- .go_cc
    testthat::expect_true(methods::is(go, "GOSubgraph"))

    old <- go@ontology
    on.exit({ go@ontology <- old }, add = TRUE)

    go@ontology <- NA_character_

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid `ontology`"
    )
})


testthat::test_that(
    ".assert_go_like() require_depth needs GO, not GOSubgraph", {
    go <- .go_cc
    testthat::expect_true(methods::is(go, "GOSubgraph"))

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go", require_depth = TRUE),
        "must be a GO object with cached depth"
    )
})


testthat::test_that(".assert_go_like() require_depth errors on empty depth", {
    testthat::skip_if_not_installed("GO.db")

    go <- GOcontext::load_go(ont = "CC", compute_depth = FALSE)
    testthat::expect_true(methods::is(go, "GO"))
    testthat::expect_true(length(go@depth) == 0L)

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go", require_depth = TRUE),
        "has no cached depth"
    )
})
