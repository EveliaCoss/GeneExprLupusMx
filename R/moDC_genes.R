#' Genes described in the differentiation of monocytes into moDCs
#'
#' Curated list of genes reported in the literature as up- or down-regulated
#' when monocytes differentiate in vitro into monocyte-derived dendritic cells
#' (moDCs), grouped by biological process. Used to check the expected
#' transcriptional changes (Figure 4C of the article).
#'
#' @format A data frame with 90 rows and 7 columns:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{Gene_name}{Name used in the source publication (e.g., \code{"IL-1beta"}).}
#'   \item{process}{Biological process (e.g., \code{"Maturation"}, \code{"Cytokines"}).}
#'   \item{Associated_pathways}{Pathway or description reported in the source.}
#'   \item{Expression}{Reported direction: \code{Upregulated} or \code{Downregulated}.}
#'   \item{link}{Link to the source publication.}
#'   \item{contrast}{Comparison in which the change was reported.}
#' }
#'
#' @details The script that generates it is \code{data-raw/lesson5_data.R}.
#'   Text was converted to ASCII (Greek letters spelled out) and the spreadsheet
#'   artifact \code{"GBP 2.00"} was corrected to \code{"GBP2"}.
#'   Eleven entries use aliases or old names (\code{CD31}, \code{CD54},
#'   \code{IL8}, \code{IRF7A}, \code{ERF2}, \code{CD32}, \code{TNFR2},
#'   \code{IP3KB}, \code{HSP27}, \code{GRP78}) or genes absent from the
#'   data (\code{CXCL8}); they do not match the gene symbols of
#'   \code{filtered_txi} and, as in the original analysis, are not shown.
#'
#' @source See the \code{link} column.
#'
#' @examples
#' table(moDC_genes$process)
"moDC_genes"
