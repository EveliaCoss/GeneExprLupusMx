#' Gene biotypes
#'
#' Biotype of each gene in \code{filtered_txi}, used to annotate differentially
#' expressed genes as protein-coding genes (PCG) or long non-coding RNAs (lncRNA).
#'
#' @format A data frame with 21,699 rows and 3 columns:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{gene_biotype}{Ensembl biotype (e.g., \code{"protein_coding"}, \code{"lncRNA"}).}
#'   \item{Group_biotype}{Grouped biotype: \code{PCG}, \code{lncRNA},
#'     \code{Pseudogene}, \code{small ncRNAs} or \code{TEC}.}
#' }
#'
#' @details The script that generates it is \code{data-raw/lesson4_annotations.R}.
#'
#' @source Ensembl gene annotation.
#'
#' @examples
#' table(biotype_summary$Group_biotype)
"biotype_summary"
