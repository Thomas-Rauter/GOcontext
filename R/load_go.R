#' Load a GO sub-ontology graph from \pkg{GO.db}
#'
#' @description
#' Loads a single Gene Ontology (GO) sub-ontology (BP, MF, or CC) from
#' the Bioconductor annotation package \pkg{GO.db} and constructs an
#' in-memory directed acyclic graph (DAG) representation.
#'
#' The returned \code{GO} object contains GO term metadata together with
#' the parent–child relationships required for graph traversal and
#' downstream restriction operations. The ontology snapshot corresponds
#' to the installed version of \pkg{GO.db}, which is recorded in the
#' returned object for reproducibility.
#'
#' @param ont \code{character(1)} GO sub-ontology to load. One of
#'   \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#'
#' @param include_obsolete \code{logical(1)} Whether obsolete GO terms
#'   should be included. If \code{FALSE} (default), obsolete terms are
#'   removed from the graph.
#'
#' @return
#' A \code{GO} S4 object representing the requested GO sub-ontology.
#' The object contains:
#' \itemize{
#'   \item GO term metadata
#'   \item parent–child relationships defining the ontology DAG
#'   \item adjacency lists for efficient graph traversal
#'   \item an empty GO-to-gene mapping table to be populated with
#'         \code{attach_org()}
#'   \item the \pkg{GO.db} version used to construct the graph
#' }
#'
#' @details
#' The ontology graph loaded by \code{load_go()} serves as the starting
#' point for GO restriction workflows implemented in this package.
#' Organism-specific mappings can be attached using
#' \code{attach_org()}, and the graph can subsequently be restricted
#' using functions such as \code{filter_mapped()} or \code{subset_go()}.
#'
#' @examples
#' go_cc <- load_go("CC")
#' class(go_cc)
#'
#' @export
load_go <- function(
        ont              = c("BP", "MF", "CC"),
        include_obsolete = FALSE
) {
    ont <- match.arg(ont)
    if (!is.logical(include_obsolete)
        || length(include_obsolete) != 1L
        || is.na(include_obsolete)) {
        rlang::abort(
            "`include_obsolete` must be a logical(1): TRUE or FALSE.",
            arg = "include_obsolete"
        )
    }

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
    adj <- .make_adj_lists(
        keep_ids = keep_ids,
        edges    = edges
        )

    methods::new(
        "GO",
        ontology = ont,
        version = as.character(utils::packageVersion("GO.db")),
        ids = term_df$go_id,
        terms = term_df,
        edges = edges,
        parents = adj$parents,
        children = adj$children,
        map = data.frame(
            go_id = character(0),
            gene_id = character(0),
            stringsAsFactors = FALSE
        )
    )
}


# Level 1 function definitions -------------------------------------------------


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
        CC = "GOCCPARENTS"
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
