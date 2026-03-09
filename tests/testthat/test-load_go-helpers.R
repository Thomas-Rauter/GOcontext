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
        base::names(edges),
        c("child", "parent")
    )
    testthat::expect_identical(base::nrow(edges), 0L)
})

testthat::test_that(".get_go_edges() drops edges when parents not in keep", {
    testthat::skip_if_not_installed("GO.db")

    parent_map <- GOcontext:::.get_parent_map("CC")

    kids <- base::names(parent_map)
    kids <- kids[!base::is.na(kids) & nzchar(kids)]
    testthat::skip_if(base::length(kids) == 0L)

    child <- kids[[1]]
    parents <- parent_map[[child]]
    parents <- base::as.character(parents)
    parents <- parents[!base::is.na(parents) & nzchar(parents)]
    testthat::skip_if(base::length(parents) == 0L)

    keep_ids <- child

    edges <- GOcontext:::.get_go_edges(
        parent_map = parent_map,
        keep_ids = keep_ids
    )

    testthat::expect_s3_class(edges, "data.frame")
    testthat::expect_identical(
        base::names(edges),
        c("child", "parent")
    )
    testthat::expect_identical(base::nrow(edges), 0L)
})

testthat::test_that(".get_go_edges() returns some edges for real keep_ids", {
    testthat::skip_if_not_installed("GO.db")

    parent_map <- GOcontext:::.get_parent_map("CC")

    kids <- base::names(parent_map)
    kids <- kids[!base::is.na(kids) & nzchar(kids)]
    testthat::skip_if(base::length(kids) == 0L)

    child <- NULL
    parent <- NULL

    for (k in kids) {
        ps <- parent_map[[k]]
        ps <- base::as.character(ps)
        ps <- ps[!base::is.na(ps) & nzchar(ps)]
        if (base::length(ps) > 0L) {
            child <- k
            parent <- ps[[1]]
            break
        }
    }

    testthat::skip_if(base::is.null(child) || base::is.null(parent))

    keep_ids <- base::unique(c(child, parent))

    edges <- GOcontext:::.get_go_edges(
        parent_map = parent_map,
        keep_ids = keep_ids
    )

    testthat::expect_s3_class(edges, "data.frame")
    testthat::expect_identical(
        base::names(edges),
        c("child", "parent")
    )

    testthat::expect_true(base::nrow(edges) >= 1L)
    testthat::expect_true(base::all(edges$child %in% keep_ids))
    testthat::expect_true(base::all(edges$parent %in% keep_ids))
})
