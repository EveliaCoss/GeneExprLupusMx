#' Z-score matrix checkpoint from the RNA-seq lesson
#'
#' Gene-wise z-scores of the variance-stabilized expression (VST) computed in
#' the lesson "De las cuentas de RNA-seq a la exploración de datos"
#' (\href{https://eveliacoss.github.io/GeneExprLupusMx/articles/intro_RNAseq_es.html}{intro_RNAseq_es}). Useful for heatmaps and clustering.
#'
#' This object is not installed with the package because of its size: get
#' it with \code{zscore_mat <- download_data("zscore_mat")} (see \code{\link{download_data}}).
#'
#' Steps: \code{vst(dds, blind = FALSE)}, removal of genes with zero variance,
#' and \code{t(scale(t(...)))} so that each gene has mean 0 and standard
#' deviation 1 across samples.
#'
#' @format A numeric matrix with 13,596 genes (rows) and 99 samples (columns).
#'
#' @details The script that generates it is \code{data-raw/lesson2_checkpoints.R}.
#'
#' @source Generated from \code{\link{dds}}.
#'
#' @seealso \code{\link{dds}}
#'
#' @examplesIf interactive()
#' zscore_mat <- download_data("zscore_mat")
#' dim(zscore_mat)
#' zscore_mat[1:4, 1:4]
#' @name zscore_mat
#' @docType data
#' @keywords datasets
#' @usage NULL
NULL
