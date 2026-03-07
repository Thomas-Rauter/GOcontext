#' Load a GO sub-ontology graph from \pkg{GO.db}
#'
#' @description
#' Constructs a reproducible, in-memory representation of a single GO
#' sub-ontology (BP, MF, or CC) using the Bioconductor annotation package
#' \pkg{GO.db}. The returned \code{GO} object contains GO term metadata and the
#' directed acyclic graph (DAG) structure required for seed expansion and
#' subgraph restriction.
#'
#' The ontology snapshot is defined by the installed version of \pkg{GO.db}.
#' The \pkg{GO.db} package version is recorded in the returned object for
#' reproducibility.
#'
#' @param ont \code{character(1)} GO sub-ontology to load. One of
#'   \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#'
#' @param include_obsolete \code{logical(1)} If \code{TRUE}, include obsolete GO
#'   terms in the returned object. If \code{FALSE} (default), obsolete terms
#'   are excluded.
#'
#' @param compute_depth \code{logical(1)} If \code{TRUE} (default), compute and
#'   store the depth of each GO term in the returned object. Precomputing depth
#'   speeds up downstream operations that rely on hierarchical distance. If
#'   \code{FALSE}, the depth vector is left empty and may be computed later
#'   when required.
#'
#' @return
#' A \code{GO} S4 object representing the requested GO sub-ontology, including
#' term metadata, DAG edges/adjacency, and the \pkg{GO.db} version used.
#'
#' @details
#' Conceptual steps performed by \code{load_go()}:
#' \enumerate{
#'   \item Validate inputs and optionally check the installed \pkg{GO.db}
#'     version.
#'   \item Extract GO term identifiers and term names from \pkg{GO.db} and
#'     retain only those belonging to \code{ont}.
#'   \item Optionally remove obsolete terms (using the obsolete mapping when
#'     available, otherwise by detecting \code{"(obsolete)"} suffixes in term
#'     names).
#'   \item Extract ontology edges for \code{ont} (parent-child relationships)
#'     from \pkg{GO.db} and keep only edges where both endpoints are retained.
#'   \item Build adjacency representations (parents/children) for DAG
#'     traversal.
#'   \item Return a \code{GO} object with embedded provenance, including the
#'     \pkg{GO.db} version.
#' }
#'
#' @examples
#' go_cc <- load_go(ont = "CC", compute_depth = FALSE)
#' class(go_cc)
#' head(go_terms(go_cc)[, c("go_id", "term")])
#'
#' @export
load_go <- function(
        ont              = c("BP", "MF", "CC"),
        include_obsolete = FALSE,
        compute_depth    = TRUE
) {
    ont <- match.arg(ont)
    .load_go_validate_inputs(
        include_obsolete = include_obsolete,
        compute_depth    = compute_depth
        )

    parent_map <- .get_parent_map(ont)
    term_df <- .get_go_term_df(
        parent_map       = parent_map,
        ont              = ont,
        include_obsolete = include_obsolete
        )
    keep_ids <- term_df$go_id

    edges <- .get_go_edges(
        parent_map = parent_map,
        keep_ids   = keep_ids
        )
    edges <- edges[!duplicated(edges), , drop = FALSE]
    adj <- .make_adj_lists(
        keep_ids = keep_ids,
        edges    = edges
        )

    depth <- integer(0)
    if (isTRUE(compute_depth)) {
        roots <- .get_roots(parents_list = adj$parents)
        depth <- .compute_depth(
            keep_ids      = keep_ids,
            children_list = adj$children,
            roots         = roots
        )
        storage.mode(depth) <- "integer"
    }

    ids <- term_df$go_id
    methods::new(
        "GO",
        ontology = ont,
        version = as.character(utils::packageVersion("GO.db")),
        ids = term_df$go_id,
        terms = term_df,
        edges = edges,
        parents = adj$parents,
        children = adj$children,
        depth = depth,
        map = data.frame(
            go_id = character(0),
            gene_id = character(0),
            stringsAsFactors = FALSE
        )
    )
}


# Level 1 function definitions -------------------------------------------------


#' Validate inputs for load_go()
#'
#' @param include_obsolete logical(1). Whether to include obsolete terms.
#' @param compute_depth logical(1). Whether to compute and store the depth
#'   vector in the returned GO object.
#'
#' @return Invisibly TRUE on success.
#' @noRd
.load_go_validate_inputs <- function(
        include_obsolete,
        compute_depth
) {
    # ---- include_obsolete
    if (!is.logical(include_obsolete) || length(include_obsolete) != 1L) {
        rlang::abort(
            message = "`include_obsolete` must be a logical(1): TRUE or FALSE.",
            arg = "include_obsolete"
        )
    }
    if (is.na(include_obsolete)) {
        rlang::abort(
            message = "`include_obsolete` must not be NA.",
            arg = "include_obsolete"
        )
    }

    # ---- compute_depth
    if (!is.logical(compute_depth) || length(compute_depth) != 1L) {
        rlang::abort(
            message = "`compute_depth` must be a logical(1): TRUE or FALSE.",
            arg = "compute_depth"
        )
    }
    if (is.na(compute_depth)) {
        rlang::abort(
            message = "`compute_depth` must not be NA.",
            arg = "compute_depth"
        )
    }

    invisible(TRUE)
}


#' @noRd
.get_parent_map <- function(ont) {
    env_name <- .get_parent_env_name(ont)
    parent_env <- get(env_name, envir = asNamespace("GO.db"))

    parent_map <- as.list(parent_env)
    parent_map <- .normalize_parent_map(parent_map)

    parent_map
}


#' Map GO ontology code to GO.db parent environment name
#'
#' @description
#' Returns the name of the GO.db environment containing parent relationships
#' for the specified GO sub-ontology. These environments store the DAG
#' structure used to derive parent–child term relationships.
#'
#' @param ont `character(1)` GO sub-ontology code. One of `"BP"`,
#'   `"MF"`, or `"CC"`.
#'
#' @return `character(1)` The name of the corresponding GO.db parent
#'   environment (e.g., `"GOBPPARENTS"`).
#'
#' @noRd
.get_parent_env_name <- function(ont) {
    switch(
        ont,
        BP = "GOBPPARENTS",
        MF = "GOMFPARENTS",
        CC = "GOCCPARENTS",
        rlang::abort("Unsupported ontology: {ont}.")
    )
}


#' Extract GO term identifiers from a parent map
#'
#' @description
#' Retrieves all GO term identifiers represented in a normalized GO
#' parent map. Both child identifiers (list names) and parent identifiers
#' (list values) are collected to ensure complete coverage of terms
#' represented in the ontology DAG.
#'
#' @param parent_map `list` Named list mapping GO identifiers (children)
#'   to character vectors of parent GO identifiers. Typically obtained
#'   via `.get_parent_map()`.
#'
#' @return `character()` A unique character vector of GO term identifiers
#'   appearing either as children or parents in `parent_map`. Only valid
#'   GO identifiers (matching `"^GO:[0-9]{7}$"`) are retained.
#'
#' @noRd
.get_term_ids_from_parents <- function(parent_map) {
    children <- names(parent_map)
    parents <- unique(unlist(parent_map, use.names = FALSE))

    ids <- unique(c(children, parents))
    ids <- ids[.is_go_id(ids)]
    ids
}


#' Retrieve obsolete GO term identifiers
#'
#' @description
#' Extracts GO term identifiers marked as obsolete in the installed
#' \pkg{GO.db} annotation package. If the obsolete mapping is unavailable
#' or cannot be accessed, an empty character vector is returned.
#'
#' @return `character()` A vector of GO IDs flagged as obsolete.
#'
#' @details
#' Obsolete terms are identified using the `GOOBSOLETE` object in
#' \pkg{GO.db}. If the object does not exist (depending on GO.db version),
#' the function safely returns `character(0)`.
#'
#' @noRd
.get_obsolete_goids <- function() {
    ns <- asNamespace("GO.db")
    if (!base::exists("GOOBSOLETE", where = ns, inherits = FALSE)) {
        return(character(0))
    }
    obs_obj <- base::get("GOOBSOLETE", envir = ns, inherits = FALSE)
    tryCatch(
        as.character(AnnotationDbi::keys(obs_obj)),
        error = function(e) character(0)
    )
}


#' Construct GO term metadata table from a parent map
#'
#' @description
#' Builds a data frame of GO term metadata for a specified GO
#' sub-ontology (BP, MF, or CC) using a precomputed parent map.
#' Term identifiers are derived from the names and values of
#' `parent_map`, and metadata are extracted from \pkg{GO.db}.
#'
#' @param parent_map `list` Named list mapping GO identifiers (children)
#'   to character vectors of parent GO identifiers. Typically obtained
#'   via `.get_parent_map()`.
#'
#' @param ont `character(1)` GO sub-ontology code. One of `"BP"`,
#'   `"MF"`, or `"CC"`. Used to populate the `ontology` column.
#'
#' @param include_obsolete `logical(1)` If `TRUE`, include obsolete
#'   GO terms. If `FALSE`, obsolete terms are removed.
#'
#' @return `data.frame` A data frame with columns:
#'   \itemize{
#'     \item `go_id` GO identifier.
#'     \item `term` GO term name.
#'     \item `ontology` Ontology namespace.
#'     \item `obsolete` Logical flag indicating obsolete status.
#'   }
#'
#' @details
#' GO identifiers are derived from both the child names and parent
#' values present in `parent_map`. Obsolete terms are identified either
#' via `GOOBSOLETE` or by detecting the suffix "(obsolete)" in term names.
#'
#' @noRd
.get_go_term_df <- function(
        parent_map,
        ont,
        include_obsolete
) {
    keep_ids <- .get_term_ids_from_parents(parent_map)

    term_vec <- AnnotationDbi::mapIds(
        x         = GO.db::GO.db,
        keys      = keep_ids,
        column    = "TERM",
        keytype   = "GOID",
        multiVals = "first",
        select    = FALSE
    )

    df <- data.frame(
        go_id            = names(term_vec),
        term             = unname(term_vec),
        ontology         = ont,
        stringsAsFactors = FALSE
    )

    obsolete_ids <- .get_obsolete_goids()
    if (length(obsolete_ids) > 0L) {
        df$obsolete <- df$go_id %in% obsolete_ids
    } else {
        df$obsolete <- grepl("\\(obsolete\\)$", df$term)
    }

    if (!include_obsolete) {
        df <- df[!df$obsolete, , drop = FALSE]
    }

    df <- df[order(df$go_id), , drop = FALSE]
    rownames(df) <- df$go_id
    df
}


#' Extract GO parent–child edges from a normalized parent map
#'
#' @description
#' Constructs directed parent–child relationships for a specified set of
#' GO identifiers from a precomputed GO parent map. The parent map is
#' typically derived from the corresponding `*PARENTS` environment in
#' \pkg{GO.db} and normalized to contain character vectors without `NA`.
#'
#' @param parent_map `list` Named list mapping GO identifiers (children)
#'   to character vectors of parent GO identifiers. Usually obtained via
#'   `.get_parent_map()`.
#'
#' @param keep_ids `character()` Vector of GO identifiers to retain in the
#'   returned edge set.
#'
#' @return `data.frame` A two-column data frame with columns:
#'   \itemize{
#'     \item `child` Child GO identifier.
#'     \item `parent` Parent GO identifier.
#'   }
#'
#' @details
#' Only edges where both `child` and `parent` are contained in
#' `keep_ids` are retained. If no such edges exist, an empty data frame
#' with the appropriate columns is returned.
#'
#' @noRd
.get_go_edges <- function(
        parent_map,
        keep_ids
) {
    .validate_keep_ids(keep_ids)

    # restrict to children we keep
    kids <- intersect(names(parent_map), keep_ids)
    if (!length(kids)) {
        return(data.frame(
            child = character(),
            parent = character(),
            stringsAsFactors = FALSE)
            )
    }

    pm <- parent_map[kids]

    pm <- lapply(pm, function(p) {
        p <- p[p %in% keep_ids]
        if (length(p)) unique(p) else character(0)
    })

    lens <- lengths(pm)
    ok <- lens > 0L
    if (!any(ok)) {
        return(data.frame(
            child = character(),
            parent = character(),
            stringsAsFactors = FALSE)
            )
    }

    pm <- pm[ok]
    kids_ok <- names(pm)

    data.frame(
        child  = rep.int(kids_ok, times = lengths(pm)),
        parent = unlist(pm, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}


#' Construct adjacency lists for GO DAG traversal
#'
#' @description
#' Builds parent and child adjacency lists from a set of GO term
#' identifiers and corresponding parent–child edges. These adjacency
#' lists facilitate efficient DAG traversal operations.
#'
#' @param keep_ids `character()` Vector of GO identifiers to include.
#'
#' @param edges `data.frame` Parent–child edge table with columns
#'   `child` and `parent`.
#'
#' @return `list` A named list with elements:
#'   \itemize{
#'     \item `parents` Named list mapping each GO ID to its parent IDs.
#'     \item `children` Named list mapping each GO ID to its child IDs.
#'   }
#'
#' @noRd
.make_adj_lists <- function(
        keep_ids,
        edges
) {
    .validate_keep_ids(keep_ids)

    # assumes edges are already deduplicated
    parents <- split(edges$parent, edges$child)
    children <- split(edges$child, edges$parent)

    # fill in missing ids with character(0) while preserving ordering
    parents <- parents[keep_ids]
    children <- children[keep_ids]

    parents <- lapply(
        parents,
        function(x) if (is.null(x)) character(0) else x
        )
    children <- lapply(
        children,
        function(x) if (is.null(x)) character(0) else x
        )

    names(parents) <- keep_ids
    names(children) <- keep_ids

    list(
        parents  = parents,
        children = children
        )
}


#' Identify root nodes in a GO subgraph
#'
#' @description
#' Determines root GO terms within a restricted set of GO identifiers.
#' A root is defined as a term that has no parents within the restricted graph.
#'
#' @param parents_list `list` Named list mapping GO IDs to their parent IDs.
#'
#' @return `character()` Vector of GO IDs representing root nodes.
#'
#' @noRd
.get_roots <- function(parents_list) {
    if (!is.list(parents_list) || is.null(names(parents_list))) {
        rlang::abort("`parents_list` must be a named list.")
    }
    names(parents_list)[lengths(parents_list) == 0L]
}


#' Computes the depth of each GO term in a directed acyclic graph
#' relative to identified root nodes.
#'
#' @description
#' Depth is defined here as the **minimum** number of edges on any path
#' from a root to the term (i.e., the shortest-path distance in an
#' unweighted DAG).
#'
#' @param keep_ids `character()` Vector of GO identifiers.
#'
#' @param children_list `list` Named list mapping GO IDs to their
#'   child identifiers.
#'
#' @param roots `character()` Vector of GO IDs designated as roots.
#'
#' @return `integer()` Named integer vector of depths, indexed by
#'   GO identifier. Root nodes have depth 0; unreachable nodes
#'   remain `NA`.
#'
#' @details
#' Depth is computed using a breadth-first traversal from all root nodes
#' simultaneously, which yields the minimum (shortest-path) distance.
#' Note that other conventions exist (e.g., maximum distance from a root);
#' this function intentionally uses the minimum-distance definition.
#'
#' @noRd
.compute_depth <- function(
        keep_ids,
        children_list,
        roots
) {
    .validate_keep_ids(keep_ids)
    depth <- stats::setNames(
        rep.int(
            NA_integer_,
            length(keep_ids)
        ),
        keep_ids
    )

    queue <- roots
    depth[roots] <- 0L

    head <- 1L
    while (head <= length(queue)) {
        current <- queue[[head]]
        head <- head + 1L

        kids <- children_list[[current]]
        if (length(kids) == 0L) next

        new_kids <- kids[is.na(depth[kids])]
        if (length(new_kids) == 0L) next

        depth[new_kids] <- depth[current] + 1L
        queue <- c(queue, new_kids)
    }

    depth
}


# Level 2 function definitions -------------------------------------------------


#' Normalize GO.db parent map values
#'
#' @param parent_map named list from GO.db *PARENTS env
#' @return named list with each element a character() vector (no NA)
#' @noRd
.normalize_parent_map <- function(parent_map) {
    parent_map <- lapply(parent_map, function(x) {
        if (is.null(x)) {
            return(character(0))
        }
        x <- as.character(unname(x))
        x <- x[!is.na(x)]
        x
    })
    parent_map
}


#' Validate GO ID vector
#'
#' @param keep_ids character() GO IDs
#' @param arg character(1) argument name for error reporting
#' @noRd
.validate_keep_ids <- function(keep_ids, arg = "keep_ids") {
    if (!is.character(keep_ids)) {
        rlang::abort(sprintf(
            "`%s` must be a character vector of GO IDs.",
            arg
            ), arg = arg)
    }
    if (anyNA(keep_ids)) {
        rlang::abort(sprintf(
            "`%s` must not contain NA values.",
            arg
            ), arg = arg)
    }
    if (anyDuplicated(keep_ids)) {
        rlang::abort(sprintf(
            "`%s` must not contain duplicated GO IDs.",
            arg
            ), arg = arg)
    }
    invisible(TRUE)
}


#' @noRd
.is_go_id <- function(x) {
    !is.na(x) & nzchar(x) & grepl("^GO:[0-9]{7}$", x)
}
