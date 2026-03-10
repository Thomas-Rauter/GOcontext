testthat::test_that(
    ".expand_to_descendants expands through children in universe",
    {
        children <- list(
            "GO:0000001" = c("GO:0000002", "GO:0000003", "GO:9999999"),
            "GO:0000002" = c("GO:0000004"),
            "GO:0000003" = character(0),
            "GO:0000004" = character(0)
        )

        universe <- c(
            "GO:0000001",
            "GO:0000002",
            "GO:0000003",
            "GO:0000004"
        )

        out <- GOcontext:::.expand_to_descendants(
            seeds = "GO:0000001",
            children = children,
            universe = universe
        )

        testthat::expect_true(is.character(out))
        testthat::expect_setequal(
            out,
            c(
                "GO:0000001",
                "GO:0000002",
                "GO:0000003",
                "GO:0000004"
            )
        )
    }
)

testthat::test_that(
    ".expand_to_descendants does not re-add already seen descendants",
    {
        children <- list(
            "GO:0000001" = c("GO:0000002", "GO:0000003"),
            "GO:0000002" = c("GO:0000003", "GO:0000004"),
            "GO:0000003" = c("GO:0000004"),
            "GO:0000004" = character(0)
        )

        universe <- c(
            "GO:0000001",
            "GO:0000002",
            "GO:0000003",
            "GO:0000004"
        )

        out <- GOcontext:::.expand_to_descendants(
            seeds = "GO:0000001",
            children = children,
            universe = universe
        )

        testthat::expect_true(is.character(out))
        testthat::expect_identical(anyDuplicated(out), 0L)
        testthat::expect_setequal(
            out,
            c(
                "GO:0000001",
                "GO:0000002",
                "GO:0000003",
                "GO:0000004"
            )
        )
    }
)
