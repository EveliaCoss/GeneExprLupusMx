#' Transcriptional signature of tolerogenic dendritic cells (tolDCs)
#'
#' Genes reported as activated (up-regulated) or repressed (down-regulated) in
#' tolerogenic dendritic cells by Robertson et al. (2021). Used to evaluate
#' whether the tolDCs of the study acquire the expected tolerogenic signature
#' (Figure 6 of the article; \code{vignette("l6_DEGs_tolDCs_es")}).
#'
#' @format A data frame with 86 rows and 4 columns:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{tolDC_signature}{Category reported by Robertson et al.: direction
#'     (\code{Activated} or \code{Repressed}) and set in which the gene was
#'     found: tolDCs generated with IL-10 (\code{IL10}), with IL-10 plus
#'     dexamethasone (\code{IL10+DXM}) or the consensus tolDC signature
#'     (\code{tolDC_signature}). Genes found in several sets list all of them,
#'     separated by commas.}
#'   \item{biotype}{Gene biotype (\code{PCG} or \code{lncRNA}).}
#'   \item{link}{Link to the source publication.}
#' }
#'
#' @details The script that generates it is \code{data-raw/lesson6_annotations.R}.
#'
#' @source Robertson H, et al. (2021). Transcriptomic analysis identifies a
#'   tolerogenic dendritic cell signature. Front Immunol 12:733231.
#'   \doi{10.3389/fimmu.2021.733231}
#'
#' @examples
#' table(tolDC_genes$tolDC_signature)
"tolDC_genes"
