#' Interferon-stimulated genes (ISGs)
#'
#' List of 423 interferon-stimulated genes, used to identify differentially
#' expressed ISGs.
#'
#' @format A data frame with 423 rows and 1 column:
#' \describe{
#'   \item{genes}{Gene symbol.}
#' }
#'
#' @details The script that generates it is \code{data-raw/lesson4_annotations.R}.
#'
#' @source Yu et al. (2025), Table S3: "The list of interferon-stimulated genes"
#'   (\doi{10.6084/m9.figshare.28353289}).
#'
#' @examples
#' head(isg_genes)
#' "IFI27" %in% isg_genes$genes
"isg_genes"
