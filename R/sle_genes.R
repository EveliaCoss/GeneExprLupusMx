#' Genes associated with systemic lupus erythematosus (SLE)
#'
#' Genes reported as associated with SLE in the literature, used to identify
#' SLE-related genes among the differentially expressed genes.
#'
#' @format A data frame with 219 rows (213 unique genes; some genes are
#'   reported by more than one study) and 4 columns:
#' \describe{
#'   \item{Gene_name}{Gene symbol.}
#'   \item{Associated_pathways}{Pathway or function associated with the gene
#'     (\code{NA} if not reported).}
#'   \item{Reference}{Study that reports the association.}
#'   \item{Link}{Link to the study.}
#' }
#'
#' @details Compiled from Shen et al. (2022), Deng and Tsao (2017),
#'   Bing et al. (2016) and Ortiz-Fernandez et al. (2023). The script that
#'   generates it is \code{data-raw/lesson4_annotations.R}.
#'
#' @source See the \code{Link} column.
#'
#' @examples
#' table(sle_genes$Reference)
#' length(unique(sle_genes$Gene_name))
"sle_genes"
