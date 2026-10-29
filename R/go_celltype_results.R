#' Functional enrichment checkpoint: dendritic cells vs monocytes (g:Profiler)
#'
#' Functional enrichment of the differentially expressed genes (DEG) between
#' dendritic cells and monocytes (\code{\link{de_celltype_results}}), obtained
#' with \pkg{gprofiler2} in the lesson about monocyte differentiation
#' (\code{vignette("l5_DEGs_monocytes_diff_es")}).
#'
#' DEG were defined with \code{|log2FoldChange| >= 1} and \code{padj < 0.05}
#' and queried as eight lists (\code{<Ctrl|SLE>_<UP|DOWN>_<moDC|tolDC>}), one at
#' a time, with \code{gost(organism = "hsapiens", correction_method = "fdr",
#' user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)}.
#'
#' @format A data frame with one row per significant term and the same 10
#'   columns as \code{\link{go_results}}: \code{query}, \code{source},
#'   \code{term_id}, \code{term_name}, \code{term_size}, \code{query_size},
#'   \code{intersection_size}, \code{p_value}, \code{parents} and
#'   \code{intersection}.
#'
#' @details The archived g:Profiler version \code{e113_eg59_p19} (2025) was
#'   used; it is stored in \code{attr(go_celltype_results, "gprofiler_version")}.
#'   The script that generates it is \code{data-raw/lesson5_checkpoints.R}.
#'
#' @source g:Profiler (\url{https://biit.cs.ut.ee/gprofiler/}).
#'
#' @seealso \code{\link{de_celltype_results}}, \code{\link{go_results}}
#'
#' @examples
#' table(go_celltype_results$query)
"go_celltype_results"
