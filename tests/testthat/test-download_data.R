# Pruebas de download_data(), sin internet: en lugar del release de GitHub se
# "descarga" un archivo local (file://), simulando release_url()

# file:// URL de un .rds local con `x`
local_release <- function(x, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".rds", .local_envir = env)
  saveRDS(x, f)
  path <- normalizePath(f, winslash = "/")
  paste0("file://", if (.Platform$OS.type == "windows") "/", path)
}

cache_file <- function(name, tag) {
  file.path(tools::R_user_dir("GeneExprLupusMx", "cache"), tag, paste0(name, ".rds"))
}

test_that("download_data only accepts the three large objects", {
  expect_error(download_data("metadata_cleaned"), "should be one of")
})

test_that("download_data downloads once and then reads the saved copy", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  fake <- matrix(1:4, 2)
  url <- local_release(fake)
  local_mocked_bindings(release_url = function(name, tag) url)

  expect_message(x <- download_data("zscore_mat", tag = "test-tag"), "Downloading zscore_mat")
  expect_identical(x, fake)
  expect_true(file.exists(cache_file("zscore_mat", "test-tag")))
  # la segunda vez no descarga: lee la copia guardada
  expect_no_message(y <- download_data("zscore_mat", tag = "test-tag"))
  expect_identical(y, fake)
})

test_that("a failed download gives a clear error and is not saved", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  local_mocked_bindings(release_url = function(name, tag) "file:///no/such/file.rds")
  expect_error(suppressMessages(download_data("dds", tag = "test-tag")),
               "Could not download 'dds'")
  expect_false(file.exists(cache_file("dds", "test-tag")))
})

test_that("download_data still works if the cache folder is not writable", {
  skip_on_os("windows") # Sys.chmod() no quita permisos de escritura en Windows
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  dir <- dirname(cache_file("dds", "test-tag"))
  dir.create(dir, recursive = TRUE)
  Sys.chmod(dir, "0555")
  withr::defer(Sys.chmod(dir, "0755"))
  skip_if(file.access(dir, 2) == 0, "the cache folder is still writable (root user?)")

  fake <- list(a = 1)
  url <- local_release(fake)
  local_mocked_bindings(release_url = function(name, tag) url)
  expect_warning(x <- suppressMessages(download_data("dds", tag = "test-tag")),
                 "Could not save 'dds' in the cache folder")
  expect_identical(x, fake)
  expect_false(file.exists(cache_file("dds", "test-tag")))
})

test_that("release_url points to the GitHub release of the package", {
  expect_identical(
    release_url("dds", "v0.1.0"),
    "https://github.com/EveliaCoss/GeneExprLupusMx/releases/download/v0.1.0/dds.rds"
  )
})
