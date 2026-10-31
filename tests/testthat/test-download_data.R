# Pruebas de download_data()

test_that("download_data only accepts the three large objects", {
  expect_error(download_data("metadata_cleaned"), "should be one of")
})

test_that("download_data reads the saved copy without downloading again", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  path <- file.path(tools::R_user_dir("GeneExprLupusMx", "cache"), "test-tag",
                    "zscore_mat.rds")
  dir.create(dirname(path), recursive = TRUE)
  fake <- matrix(1:4, 2)
  saveRDS(fake, path)
  expect_identical(download_data("zscore_mat", tag = "test-tag"), fake)
})

test_that("a failed download gives a clear error and is not saved", {
  skip_on_cran()
  skip_if_offline()
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  expect_error(suppressMessages(download_data("dds", tag = "no-such-release")),
               "Could not download 'dds'")
  path <- file.path(tools::R_user_dir("GeneExprLupusMx", "cache"),
                    "no-such-release", "dds.rds")
  expect_false(file.exists(path))
})
