#' Example GO cellular component graph
#'
#' A Gene Ontology (GO) graph restricted to the Cellular Component (CC)
#' ontology. This dataset is intended for use in examples and documentation.
#'
#' The object is of class \code{GO} and contains GO term metadata and
#' parent–child relationships for the CC ontology.
#'
#' @format An S4 \code{GO} object containing term metadata and the
#' parent–child DAG structure.
#'
#' @details
#' The graph was created with:
#'
#' \preformatted{
#' load_go(ont = "CC")
#' }
#'
#' @source Gene Ontology
#'
#' @examples
#' data("go_cc", package = "GOcontext")
#' go_cc
#'
#' @docType data
#' @name go_cc
#' @keywords internal
NULL


#' Example GO CC graph with E. coli gene mappings
#'
#' A Gene Ontology Cellular Component graph with organism mappings
#' attached for \emph{Escherichia coli} K-12 using the
#' \code{org.EcK12.eg.db} annotation package.
#'
#' This dataset is intended for use in runnable examples demonstrating
#' functions such as \code{\link{attach_org}}, \code{\link{filter_mapped}},
#' \code{\link{as_term2gene}}, and \code{\link{as_term2name}}.
#'
#' @format An S4 \code{GO} object with an attached GO-to-gene mapping.
#'
#' @details
#' The object was generated with:
#'
#' \preformatted{
#' go_cc <- load_go(ont = "CC")
#' go_cc_ecoli <- attach_org(
#'     go = go_cc,
#'     OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
#'     keytype = "ENTREZID"
#' )
#' }
#'
#' @source Gene Ontology and org.EcK12.eg.db
#'
#' @examples
#' data("go_cc_ecoli", package = "GOcontext")
#' go_cc_ecoli
#'
#' @docType data
#' @name go_cc_ecoli
#' @keywords internal
NULL
