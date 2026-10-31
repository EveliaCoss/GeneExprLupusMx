# Pruebas de get_heatmap_mat() y create_heatmap()

# DESeqDataSet pequeño y simulado: 1,200 genes (vst() necesita al menos 1,000)
# y un gen sin cuentas, que tiene varianza cero y debe eliminarse
make_small_dds <- function() {
  set.seed(1)
  counts <- matrix(stats::rnbinom(1200 * 6, mu = 100, size = 10), ncol = 6,
                   dimnames = list(paste0("gene", 1:1200), paste0("S", 1:6)))
  counts["gene1", ] <- 0L
  col_data <- data.frame(Group = factor(rep(c("Ctrl", "SLE"), each = 3)),
                         row.names = colnames(counts))
  suppressMessages(
    DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ Group)
  )
}

# Matriz pequeña de Z-scores simulados: 50 genes x 20 muestras
make_small_zscores <- function() {
  set.seed(1)
  mat <- matrix(stats::rnorm(50 * 20), nrow = 50,
                dimnames = list(paste0("gene", 1:50), metadata_cleaned$sample_ID[1:20]))
  t(scale(t(mat)))
}

test_that("get_heatmap_mat returns Z-scores and drops zero-variance genes", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  small_dds <- make_small_dds()
  mat <- suppressMessages(get_heatmap_mat(small_dds))

  expect_true(is.matrix(mat))
  expect_identical(colnames(mat), colnames(small_dds))
  expect_false("gene1" %in% rownames(mat))
  expect_identical(nrow(mat), 1199L)
  expect_lt(max(abs(rowMeans(mat))), 1e-8)
  expect_equal(unname(apply(mat, 1, stats::sd)), rep(1, nrow(mat)), tolerance = 1e-8)
})

test_that("get_heatmap_mat orders the columns with sample_order", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  small_dds <- make_small_dds()
  new_order <- rev(colnames(small_dds))
  mat <- suppressMessages(get_heatmap_mat(small_dds, sample_order = new_order))
  expect_identical(colnames(mat), new_order)
})

test_that("get_heatmap_mat works on one cell type of dds", {
  skip_on_cran()
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  dds <- local_download("dds")
  mo <- dds[, dds$Cell_type == "monocyte"]
  DESeq2::design(mo) <- ~ Group
  mat <- suppressMessages(get_heatmap_mat(mo))
  expect_identical(ncol(mat), 33L)
  expect_true(all(rownames(mat) %in% rownames(dds)))
})

test_that("create_heatmap returns a ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  mat <- make_small_zscores()
  ht <- create_heatmap(mat, color = "black", split = NULL)

  expect_s4_class(ht, "Heatmap")
  expect_identical(dim(ht@matrix), dim(mat))
  expect_identical(ht@name, "Normalized gene\nexpression (Z-score)")
  expect_s4_class(ht@right_annotation, "HeatmapAnnotation")
})

test_that("create_heatmap splits columns and can be drawn and stacked", {
  skip_if_not_installed("ComplexHeatmap")
  mat <- make_small_zscores()
  group <- metadata_cleaned$Group[1:20]
  ht_mo <- create_heatmap(mat, color = "black", split = group)
  ht_dc <- create_heatmap(mat, color = "palevioletred", split = group)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(ComplexHeatmap::draw(ht_mo))
  # todos los heatmaps comparten el mismo nombre de leyenda; ComplexHeatmap avisa
  expect_no_error(suppressWarnings(ComplexHeatmap::draw(ComplexHeatmap::`%v%`(ht_mo, ht_dc))))
})
