#' Create a ComplexHeatmap with row annotation
#'
#' Generates a heatmap using the \pkg{ComplexHeatmap} package with a single-color
#' row annotation. Designed to stack several heatmaps (for example, one per cell
#' type) with consistent formatting.
#'
#' @param mat A numeric matrix representing the heatmap values (e.g., Z-scores).
#' @param color A character string specifying the annotation color
#'   (e.g., \code{"black"}, \code{"gold1"}).
#' @param split A factor, character vector or data frame used to split the
#'   columns of the heatmap. \code{NULL} for no split.
#' @param common_names A character vector of row names, used to size the
#'   annotation. By default, the row names of \code{mat}.
#' @param use_raster_flag Logical. If \code{TRUE}, rasterization is enabled for
#'   faster rendering of large matrices (requires the \pkg{magick} package and a
#'   working Cairo/X11 installation). Default \code{FALSE}.
#'
#' @return A \code{Heatmap} object from \pkg{ComplexHeatmap}, ready to be
#'   combined or drawn.
#'
#' @details Requires the Bioconductor package \pkg{ComplexHeatmap}.
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#'   # small simulated Z-score matrix: 50 genes x 20 samples
#'   # (with the lesson data: zscore_mat <- download_data("zscore_mat"))
#'   set.seed(1)
#'   mat <- t(scale(t(matrix(rnorm(50 * 20), nrow = 50))))
#'   ht <- create_heatmap(mat, color = "black", split = NULL)
#'   ComplexHeatmap::draw(ht)
#' }
#'
#' @export
create_heatmap <- function(mat, color, split, common_names = rownames(mat),
                           use_raster_flag = FALSE) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Install ComplexHeatmap to use this function: BiocManager::install('ComplexHeatmap')")
  }
  ha <- ComplexHeatmap::rowAnnotation(
    m = rep("1", length(common_names)),
    show_legend = FALSE,
    col = list(m = c("1" = color)),
    show_annotation_name = FALSE,
    annotation_label = NULL
  )

  ComplexHeatmap::Heatmap(
    mat,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = TRUE,
    name = "Normalized gene\nexpression (Z-score)",
    heatmap_legend_param = list(direction = "horizontal"),
    column_split = split,
    right_annotation = ha,
    use_raster = use_raster_flag
  )
}
