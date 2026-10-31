#' Download the large data objects of the package
#'
#' The three largest objects (\code{\link{filtered_txi}}, \code{\link{dds}} and
#' \code{\link{zscore_mat}}) are not installed with the package, to keep it
#' small. They are stored as files of the GitHub release of the package, and
#' this function downloads them the first time they are used. Later calls read
#' the copy saved on your computer, so no internet connection is needed.
#'
#' @param name Name of the object: \code{"filtered_txi"}, \code{"dds"} or
#'   \code{"zscore_mat"}.
#' @param tag GitHub release that contains the files. The default is the
#'   release with the data used in the lessons.
#'
#' @return The requested object.
#'
#' @details Files are saved in \code{tools::R_user_dir("GeneExprLupusMx",
#'   "cache")}. Delete that folder to free space or to download them again.
#'   The first download needs about 46 MB in total.
#'
#' @seealso \code{\link{filtered_txi}}, \code{\link{dds}},
#'   \code{\link{zscore_mat}}
#'
#' @examplesIf interactive()
#' filtered_txi <- download_data("filtered_txi")
#' dim(filtered_txi$counts)
#'
#' @export
download_data <- function(name = c("filtered_txi", "dds", "zscore_mat"),
                          tag = "v0.1.0") {
  name <- match.arg(name)
  path <- file.path(tools::R_user_dir("GeneExprLupusMx", "cache"), tag,
                    paste0(name, ".rds"))
  if (!file.exists(path)) {
    url <- paste0("https://github.com/EveliaCoss/GeneExprLupusMx/releases/download/",
                  tag, "/", name, ".rds")
    message("Downloading ", name, " (only the first time)...")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    old <- options(timeout = max(600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    ok <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      stop("Could not download '", name, "' from ", url,
           ". Check your internet connection and try again.", call. = FALSE)
    }
    # copy only complete downloads, so a failed download is never cached;
    # if the cache folder is not writable (e.g., on a server), read the download
    if (!file.copy(tmp, path, overwrite = TRUE)) return(readRDS(tmp))
  }
  readRDS(path)
}
