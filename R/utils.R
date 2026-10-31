# Packages of `pkgs` that are not installed (separate function so the tests
# can simulate a missing package)
#' @noRd
pkgs_missing <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
}

# Stops with a message that says which packages are missing and how to
# install them (BiocManager installs both Bioconductor and CRAN packages)
#' @noRd
check_pkgs <- function(pkgs, fun) {
  missing <- pkgs_missing(pkgs)
  if (length(missing) > 0) {
    stop(fun, "() needs the package(s) ", paste0("'", missing, "'", collapse = ", "),
         ", which are not installed. Install them with:\n",
         "  install.packages(\"BiocManager\")\n",
         "  BiocManager::install(c(", paste0("\"", missing, "\"", collapse = ", "), "))",
         call. = FALSE)
  }
  invisible(TRUE)
}
