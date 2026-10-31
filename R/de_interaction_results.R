#' Differential expression results: between-group comparison of contrasts
#'
#' \pkg{DESeq2} interaction results that test whether the change from monocytes
#' to in vitro differentiated dendritic cells (moDC or tolDC) is different in
#' SLE patients than in healthy controls. Used in the lesson about the
#' between-group comparison of contrasts
#' (\href{https://eveliacoss.github.io/GeneExprLupusMx/articles/l7_DEGs_between_groups_es.html}{l7_DEGs_between_groups_es}; Figure 7A and 7B of the article).
#'
#' Design (all genes of \code{filtered_txi}, the 99 samples, reference levels
#' \code{Ctrl} and \code{monocyte}): \code{~ Group + Cell_type + Group:Cell_type}.
#' Results were extracted with \code{results(alpha = 0.05, name =
#' "GroupSLE.Cell_type<moDC|tolDC>")} (Wald test, Benjamini-Hochberg adjustment).
#' The \code{log2FoldChange} is a difference of differences:
#' (dendritic cell vs monocyte in SLE) minus (dendritic cell vs monocyte in Ctrl).
#'
#' @format A named list with two data frames (\code{moDC} and \code{tolDC}),
#'   each with 21,718 rows (genes) and 7 columns: \code{gene}, \code{baseMean},
#'   \code{log2FoldChange}, \code{lfcSE}, \code{stat}, \code{pvalue} and
#'   \code{padj} (see \code{\link{de_results}}).
#'
#' @details With \code{padj < 0.05} (up / down): moDC 1 / 28 and
#'   tolDC 3 / 38. The script that generates this object is
#'   \code{data-raw/lesson7_checkpoints.R}.
#'
#' @source Generated from \code{\link{filtered_txi}} and \code{\link{metadata_cleaned}}.
#'
#' @seealso \code{\link{de_celltype_results}} (dendritic cells vs monocytes
#'   within each group), \code{\link{go_interaction_results}}.
#'
#' @examples
#' names(de_interaction_results)
#' head(de_interaction_results$tolDC[order(de_interaction_results$tolDC$padj), ])
"de_interaction_results"
