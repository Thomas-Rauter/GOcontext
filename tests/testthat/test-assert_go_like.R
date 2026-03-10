# tests/testthat/test-assert_go_like.R

testthat::test_that(".assert_go_like() accepts a valid GO object", {
    testthat::expect_invisible(
        GOcontext:::.assert_go_like(go_cc, arg = "go")
    )
})

testthat::test_that(".assert_go_like() accepts a valid GOSubgraph object", {
    go_sub <- GOcontext::subset_go(
        go = go_cc,
        ids = go_cc@terms$go_id[[1]],
        mode = "keep"
    )

    testthat::expect_true(methods::is(go_sub, "GOSubgraph"))

    testthat::expect_invisible(
        GOcontext:::.assert_go_like(go_sub, arg = "go")
    )
})

testthat::test_that(".assert_go_like() errors on invalid class", {
    testthat::expect_error(
        GOcontext:::.assert_go_like("not_a_go_object", arg = "go"),
        "must be a GO or GOSubgraph object"
    )
})

testthat::test_that(".assert_go_like() errors on invalid ontology slot", {
    go <- go_cc
    old <- go@ontology
    on.exit({ go@ontology <- old }, add = TRUE)

    go@ontology <- NA_character_

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `ontology`"
    )
})

testthat::test_that(".assert_go_like() errors on invalid version slot", {
    go <- go_cc
    old <- go@version
    on.exit({ go@version <- old }, add = TRUE)

    go@version <- NA_character_

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `version`"
    )
})

testthat::test_that(".assert_go_like() errors on invalid terms slot", {
    go <- go_cc
    old <- go@terms
    on.exit({ go@terms <- old }, add = TRUE)

    go@terms <- data.frame(
        wrong = character(0),
        stringsAsFactors = FALSE
    )

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `terms`"
    )
})

testthat::test_that(".assert_go_like() errors on invalid edges slot", {
    go <- go_cc
    old <- go@edges
    on.exit({ go@edges <- old }, add = TRUE)

    go@edges <- data.frame(
        wrong = character(0),
        stringsAsFactors = FALSE
    )

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `edges`"
    )
})

testthat::test_that(".assert_go_like() errors on invalid parents slot", {
    go <- go_cc
    old <- go@parents
    on.exit({ go@parents <- old }, add = TRUE)

    go@parents <- list(character(0))

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `parents`"
    )
})

testthat::test_that(".assert_go_like() errors on invalid children slot", {
    go <- go_cc
    old <- go@children
    on.exit({ go@children <- old }, add = TRUE)

    go@children <- list(character(0))

    testthat::expect_error(
        GOcontext:::.assert_go_like(go, arg = "go"),
        "has invalid slot `children`"
    )
})

testthat::test_that(
    ".assert_go_like() require_map errors on empty attached map",
    {
        testthat::expect_error(
            GOcontext:::.assert_go_like(
                go_cc,
                arg = "go",
                require_map = TRUE
            ),
            "must have a non-empty attached organism mapping"
        )
    }
)

testthat::test_that(
    ".assert_go_like() require_map accepts a non-empty attached map",
    {
        testthat::skip_if_not_installed("org.EcK12.eg.db")

        go_mapped <- GOcontext::attach_org(
            go = go_cc,
            OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
            keytype = "ENTREZID"
        )

        testthat::expect_invisible(
            GOcontext:::.assert_go_like(
                go_mapped,
                arg = "go",
                require_map = TRUE
            )
        )
    }
)
