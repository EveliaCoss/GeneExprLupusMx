#' Generate a heatmap-ready matrix from a DESeq2 object
#'
#' Applies the variance stabilizing transformation (VST) to a \code{DESeqDataSet},
#' removes genes with zero variance, converts each gene to Z-scores across
#' samples and orders the columns.
#'
#' @param dds_object A \code{DESeqDataSet} (for example, the samples of a
#'   single cell type).
#' @param sample_order Character vector with the sample IDs in the order in
#'   which the columns should appear. By default, the order of \code{dds_object}.
#'
#' @return A numeric matrix of Z-scores (genes in rows, samples in columns).
#'
#' @details Requires the Bioconductor packages \pkg{DESeq2} and
#'   \pkg{SummarizedExperiment}.
#'
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   mo <- dds[, dds$Cell_type == "monocyte"]
#'   DESeq2::design(mo) <- ~ Group # a single cell type: compare groups only
#'   mo_mat <- get_heatmap_mat(mo)
#'   dim(mo_mat)
#' }
#'
#' @export
get_heatmap_mat <- function(dds_object, sample_order = colnames(dds_object)) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Install DESeq2 to use this function: BiocManager::install('DESeq2')")
  }
  vsd <- DESeq2::vst(dds_object, blind = FALSE)
  mat <- SummarizedExperiment::assay(vsd)
  mat <- mat[apply(mat, 1, stats::var) != 0, ] # filter out genes with zero variance
  zscores <- t(scale(t(mat)))
  zscores[, sample_order]
}
