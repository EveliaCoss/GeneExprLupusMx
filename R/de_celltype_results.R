#' Differential expression results: dendritic cells vs monocytes
#'
#' \pkg{DESeq2} results comparing in vitro differentiated dendritic cells
#' (moDC or tolDC) against monocytes, within each group (healthy controls and
#' SLE patients). Used in the lesson about monocyte differentiation
#' (\href{https://eveliacoss.github.io/GeneExprLupusMx/articles/l5_DEGs_monocytes_diff_es.html}{l5_DEGs_monocytes_diff_es}; Figures 4 and 5 of the article).
#'
#' Designs (all genes of \code{filtered_txi}, reference level \code{monocyte}):
#' \itemize{
#'   \item Controls: \code{~ Cell_type}.
#'   \item SLE patients: \code{~ SLEDAI_category + Cell_type + Dose_category}
#'     (adjusted for disease activity and glucocorticoid dose).
#' }
#' Results were extracted with \code{results(alpha = 0.05, name =
#' "Cell_type_<moDC|tolDC>_vs_monocyte")} (Wald test, Benjamini-Hochberg adjustment).
#'
#' @format A named list with four data frames (\code{Ctrl_moDC},
#'   \code{Ctrl_tolDC}, \code{SLE_moDC}, \code{SLE_tolDC}), each with 21,718
#'   rows (genes) and 7 columns: \code{gene}, \code{baseMean},
#'   \code{log2FoldChange}, \code{lfcSE}, \code{stat}, \code{pvalue} and
#'   \code{padj} (see \code{\link{de_results}}).
#'
#' @details With \code{|log2FoldChange| >= 1} and \code{padj < 0.05}
#'   (up / down): Ctrl moDC 2688 / 2947, Ctrl tolDC 2731 / 2849,
#'   SLE moDC 2945 / 3652 and SLE tolDC 2804 / 3596. The article (DESeq2 on
#'   R 4.0.2) reports the same numbers for moDC and differs by a few genes in
#'   tolDC (2736 / 2851 and 2804 / 3592). The script that generates this object
#'   is \code{data-raw/lesson5_checkpoints.R}.
#'
#' @source Generated from \code{\link{filtered_txi}} and \code{\link{metadata_cleaned}}.
#'
#' @seealso \code{\link{de_results}} (SLE vs Ctrl within each cell type).
#'
#' @examples
#' names(de_celltype_results)
#' head(de_celltype_results$Ctrl_moDC[order(de_celltype_results$Ctrl_moDC$padj), ])
"de_celltype_results"
