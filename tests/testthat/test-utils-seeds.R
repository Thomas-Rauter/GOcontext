testthat::test_that(
    ".read_seed_catalog_rds() errors when catalog is malformed", {
    testthat::with_mocked_bindings(
        .seed_catalog_path = function(package) "fakepath",
        .seed_read_rds = function(path) list(bundles = "nope"),
        {
            testthat::expect_error(
                GOcontext:::.read_seed_catalog_rds(package = "GOcontext"),
                regexp = "malformed",
                class = "rlang_error"
            )
        },
        .package = "GOcontext"
    )
})

testthat::test_that(
    ".read_seed_catalog_rds() errors when catalog file missing", {
    testthat::with_mocked_bindings(
        .seed_catalog_path = function(package) "",
        {
            testthat::expect_error(
                GOcontext:::.read_seed_catalog_rds(package = "GOcontext"),
                regexp = "Compiled seed catalog not found",
                class = "rlang_error"
            )
        },
        .package = "GOcontext"
    )
})

testthat::test_that(".read_seed_catalog_rds() validates 'package' argument", {
    # not character(1)
    testthat::expect_error(
        GOcontext:::.read_seed_catalog_rds(package = 1),
        regexp = "package",
        class = "rlang_error"
    )

    # length != 1
    testthat::expect_error(
        GOcontext:::.read_seed_catalog_rds(package = c("GOcontext", "x")),
        regexp = "package",
        class = "rlang_error"
    )

    # NA
    testthat::expect_error(
        GOcontext:::.read_seed_catalog_rds(package = NA_character_),
        regexp = "package",
        class = "rlang_error"
    )
})

testthat::test_that(".seed_catalog_as_df() errors when bundles list is empty", {
    cfg <- list(bundles = list())  # length == 0 triggers the branch
    testthat::expect_error(
        GOcontext:::.seed_catalog_as_df(cfg),
        regexp = "No seed bundles found",
        class = "rlang_error"
    )
})

testthat::test_that(
    ".seed_bundle_to_df_row() errors when required fields missing", {
    # Missing go_ids (and others if you want)
    b <- list(
        name = "x",
        ont = "BP",
        title = "t",
        description = "d"
        # go_ids missing
    )

    testthat::expect_error(
        GOcontext:::.seed_bundle_to_df_row(b),
        regexp = "missing required field\\(s\\)",
        class = "rlang_error"
    )
})
