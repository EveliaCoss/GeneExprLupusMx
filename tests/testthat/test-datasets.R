test_that("filtered_txi elements have expected structure", {
  data("filtered_txi", package = "GeneExprLupusMx")

  # Validar que es una lista con los nombres correctos
  expect_type(filtered_txi, "list")
  expect_named(filtered_txi, c("abundance", "counts", "length", "countsFromAbundance"))

  # Validar que los elementos son matrices
  expect_true(is.matrix(filtered_txi$abundance))
  expect_true(is.matrix(filtered_txi$counts))
  expect_true(is.matrix(filtered_txi$length))

  # Validar que countsFromAbundance es un carácter
  expect_type(filtered_txi$countsFromAbundance, "character")

  # Validar dimensiones mínimas
  expect_true(nrow(filtered_txi$counts) > 0)
  expect_true(ncol(filtered_txi$counts) > 0)
})


test_that("metadata_cleaned has expected variables", {
  expect_s3_class(metadata_cleaned, "data.frame")
  expect_true("sample_ID" %in% colnames(metadata_cleaned))
  expect_true("Cell_type" %in% colnames(metadata_cleaned))
})
