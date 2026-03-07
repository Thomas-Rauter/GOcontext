# tests/testthat/test-load_go-validate_inputs.R

testthat::test_that(".load_go_validate_inputs() accepts valid inputs", {
    testthat::expect_invisible(
        GOcontext:::.load_go_validate_inputs(
            version = NULL,
            include_obsolete = FALSE,
            compute_depth = TRUE
        )
    )

    testthat::expect_invisible(
        GOcontext:::.load_go_validate_inputs(
            version = "3.19.1",
            include_obsolete = TRUE,
            compute_depth = FALSE
        )
    )

    testthat::expect_invisible(
        GOcontext:::.load_go_validate_inputs(
            version = " 3.19.1 ",
            include_obsolete = FALSE,
            compute_depth = TRUE
        )
    )
})


testthat::test_that(".load_go_validate_inputs() rejects bad version", {
    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = 1,
            include_obsolete = FALSE,
            compute_depth = TRUE
        ),
        "version.*character\\(1\\)"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = c("3.19.1", "3.19.2"),
            include_obsolete = FALSE,
            compute_depth = TRUE
        ),
        "must be length 1"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = NA_character_,
            include_obsolete = FALSE,
            compute_depth = TRUE
        ),
        "must not be NA"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = "   ",
            include_obsolete = FALSE,
            compute_depth = TRUE
        ),
        "must not be an empty string"
    )
})


testthat::test_that(".load_go_validate_inputs() rejects bad flags", {
    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = NULL,
            include_obsolete = 1,
            compute_depth = TRUE
        ),
        "include_obsolete.*logical\\(1\\)"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = NULL,
            include_obsolete = NA,
            compute_depth = TRUE
        ),
        "include_obsolete.*must not be NA"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = NULL,
            include_obsolete = FALSE,
            compute_depth = 1
        ),
        "compute_depth.*logical\\(1\\)"
    )

    testthat::expect_error(
        GOcontext:::.load_go_validate_inputs(
            version = NULL,
            include_obsolete = FALSE,
            compute_depth = NA
        ),
        "compute_depth.*must not be NA"
    )
})


# tests/testthat/test-load_go-edges-helpers.R

testthat::test_that(".validate_keep_ids() rejects invalid keep_ids", {
    testthat::expect_error(
        GOcontext:::.validate_keep_ids(keep_ids = 1),
        "must be a character vector"
    )

    testthat::expect_error(
        GOcontext:::.validate_keep_ids(
            keep_ids = c("GO:0000001", NA_character_)
        ),
        "must not contain NA"
    )

    testthat::expect_error(
        GOcontext:::.validate_keep_ids(
            keep_ids = c("GO:0000001", "GO:0000001")
        ),
        "must not contain duplicated"
    )

    testthat::expect_invisible(
        GOcontext:::.validate_keep_ids(
            keep_ids = c("GO:0000001", "GO:0000002")
        )
    )
})


testthat::test_that(".get_go_edges() returns empty if no kids kept", {
    testthat::skip_if_not_installed("GO.db")

    parent_map <- GOcontext:::.get_parent_map("CC")

    keep_ids <- c("GO:9999999", "GO:8888888")
    edges <- GOcontext:::.get_go_edges(
        parent_map = parent_map,
        keep_ids = keep_ids
        )

    testthat::expect_s3_class(edges, "data.frame")
    testthat::expect_identical(
        names(edges),
        c("child", "parent")
    )
    testthat::expect_identical(nrow(edges), 0L)
})


testthat::test_that(".get_go_edges() drops edges when parents not in keep", {
    testthat::skip_if_not_installed("GO.db")

    parent_map <- GOcontext:::.get_parent_map("CC")

    kids <- names(parent_map)
    kids <- kids[!is.na(kids) & nzchar(kids)]
    testthat::skip_if(length(kids) == 0L)

    child <- kids[[1]]
    parents <- parent_map[[child]]
    parents <- as.character(parents)
    parents <- parents[!is.na(parents) & nzchar(parents)]
    testthat::skip_if(length(parents) == 0L)

    keep_ids <- child

    edges <- GOcontext:::.get_go_edges(
        parent_map = parent_map,
        keep_ids = keep_ids
        )

    testthat::expect_s3_class(edges, "data.frame")
    testthat::expect_identical(
        names(edges),
        c("child", "parent")
    )
    testthat::expect_identical(nrow(edges), 0L)
})


testthat::test_that(".get_go_edges() returns some edges for real keep_ids", {
    testthat::skip_if_not_installed("GO.db")

    parent_map <- GOcontext:::.get_parent_map("CC")

    kids <- names(parent_map)
    kids <- kids[!is.na(kids) & nzchar(kids)]
    testthat::skip_if(length(kids) == 0L)

    child <- NULL
    parent <- NULL

    for (k in kids) {
        ps <- parent_map[[k]]
        ps <- as.character(ps)
        ps <- ps[!is.na(ps) & nzchar(ps)]
        if (length(ps) > 0L) {
            child <- k
            parent <- ps[[1]]
            break
        }
    }

    testthat::skip_if(is.null(child) || is.null(parent))

    keep_ids <- unique(c(child, parent))

    edges <- GOcontext:::.get_go_edges(
        parent_map = parent_map,
        keep_ids = keep_ids
        )

    testthat::expect_s3_class(edges, "data.frame")
    testthat::expect_identical(
        names(edges),
        c("child", "parent")
    )

    testthat::expect_true(nrow(edges) >= 1L)
    testthat::expect_true(all(edges$child %in% keep_ids))
    testthat::expect_true(all(edges$parent %in% keep_ids))
})
