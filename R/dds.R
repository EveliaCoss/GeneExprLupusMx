#' DESeqDataSet checkpoint from the RNA-seq lesson
#'
#' \code{DESeqDataSet} built in the lesson "De las cuentas de RNA-seq a la
#' exploración de datos" (\href{https://eveliacoss.github.io/GeneExprLupusMx/articles/intro_RNAseq_es.html}{intro_RNAseq_es}), ready to be used
#' as the starting point of later lessons.
#'
#' This object is not installed with the package because of its size: get
#' it with \code{dds <- download_data("dds")} (see \code{\link{download_data}}).
#'
#' Steps: \code{DESeqDataSetFromTximport()} from \code{filtered_txi} and
#' \code{metadata_cleaned} with design \code{~ Cell_type + Group} (reference
#' levels \code{monocyte} and \code{Ctrl}), removal of genes with fewer than 10
#' reads in at least 10 samples, and \code{estimateSizeFactors()}.
#'
#' @format A \code{DESeqDataSet} with 13,596 genes and 99 samples. Its
#'   \code{colData} is \code{metadata_cleaned} plus a \code{Volunteer_ID} column.
#'
#' @details Requires the \pkg{DESeq2} package (Bioconductor) to be used.
#'   The script that generates it is \code{data-raw/lesson2_checkpoints.R}.
#'
#' @source Generated from \code{filtered_txi} and \code{metadata_cleaned}.
#'
#' @seealso \code{\link{zscore_mat}}, \code{\link{filtered_txi}},
#'   \code{\link{metadata_cleaned}}
#'
#' @examplesIf interactive() && requireNamespace("DESeq2", quietly = TRUE)
#' library(DESeq2)
#' dds <- download_data("dds")
#' dds
#' design(dds)
#' @name dds
#' @docType data
#' @keywords datasets
#' @usage NULL
NULL
