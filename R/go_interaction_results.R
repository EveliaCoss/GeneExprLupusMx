#' Functional enrichment checkpoint: between-group comparison of contrasts (g:Profiler)
#'
#' Functional enrichment of the genes whose change from monocytes to dendritic
#' cells differs between SLE patients and controls
#' (\code{\link{de_interaction_results}}), obtained with \pkg{gprofiler2} in
#' \code{vignette("l7_DEGs_between_groups_es")}.
#'
#' Genes with \code{padj < 0.05} were split by the sign of
#' \code{log2FoldChange} and queried as four lists
#' (\code{SLEvsCtrl_<UP|DOWN>_<moDC|tolDC>}) with \code{gost(organism =
#' "hsapiens", correction_method = "fdr", user_threshold = 0.05, sources =
#' c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)}.
#'
#' @format A data frame with one row per significant term and the same 10
#'   columns as \code{\link{go_results}}: \code{query}, \code{source},
#'   \code{term_id}, \code{term_name}, \code{term_size}, \code{query_size},
#'   \code{intersection_size}, \code{p_value}, \code{parents} and
#'   \code{intersection}.
#'
#' @details The archived g:Profiler version \code{e113_eg59_p19} (2025) was
#'   used; it is stored in \code{attr(go_interaction_results, "gprofiler_version")}.
#'   The script that generates it is \code{data-raw/lesson7_checkpoints.R}.
#'
#' @source g:Profiler (\url{https://biit.cs.ut.ee/gprofiler/}).
#'
#' @seealso \code{\link{de_interaction_results}}, \code{\link{go_celltype_results}}
#'
#' @examples
#' table(go_interaction_results$query, go_interaction_results$source)
"go_interaction_results"
