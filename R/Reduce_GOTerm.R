# --------- processGO (function) / REDUCE GO TERMS by file----------------------
#' Process GO Enrichment Result from a Single File
#'
#' Load and analyze a single \code{.RData} file containing GO enrichment results
#' from either UP or DOWN regulated genes. Includes validation, semantic
#' similarity reduction, merging metadata, and plotting.
#'
#' @param f Character. Path to the \code{.RData} file (must include full or
#'   relative path). The file must contain exactly one data frame with the
#'   columns \code{go_id}, \code{term} and \code{p.adjust}.
#' @param direction Character. Either \code{"UP"} or \code{"DOWN"} to specify
#'   the type of enrichment.
#' @param outdir Character. Directory to save \code{.RData} output files.
#' @param figdir Character. Directory to save PNG plots.
#'
#' @details Requires the Bioconductor packages \pkg{rrvgo} and
#'   \pkg{org.Hs.eg.db}. Datasets with fewer than 40 unique GO terms are skipped.
#'
#' @return Invisibly, \code{NULL}. Saves processed files and figures in
#'   \code{outdir} and \code{figdir}.
#'
#' @seealso \code{\link{Reduce_GOTerm}} for batch processing of multiple
#'   enrichment files.
#'
#' @examples
#' \dontrun{
#' # Analyze a specific UP-regulated dataset
#' processGO(
#'   f = "data/enrichment/UP_GO_MySample1.RData",
#'   direction = "UP",
#'   outdir = "results/GO_terms/",
#'   figdir = "figures/GO_terms/"
#' )
#' }
#'
#' @export
processGO <- function(f, direction = c("UP", "DOWN"), outdir, figdir) {
  if (!requireNamespace("rrvgo", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install rrvgo and org.Hs.eg.db to use this function: ",
         "BiocManager::install(c('rrvgo', 'org.Hs.eg.db'))")
  }
  direction <- match.arg(direction)
  plot_name <- gsub(paste0("^", direction, "_GO_(.+)\\.RData$"), "\\1", basename(f))

  message("\n------------------------------------")
  message("Dataset: ", plot_name)
  message("------------------------------------")

  # Cargar el archivo y detectar el nombre del objeto
  env <- new.env() # entorno aparte, sin contaminar el resto
  load(f, envir = env)
  obj_names <- ls(env)

  if (length(obj_names) != 1) {
    stop("El archivo debe contener exactamente un objeto. Se encontraron: ",
         paste(obj_names, collapse = ", "))
  }

  bar_data <- env[[obj_names[1]]]

  if (!("go_id" %in% colnames(bar_data))) {
    message("No 'go_id' column found. Skipping file: ", f)
    return(invisible(NULL))
  }

  bar_data <- bar_data[!is.na(bar_data$go_id), ]
  bar_data <- bar_data[!duplicated(bar_data$go_id), ]

  if (nrow(bar_data) < 40) {
    message("Dataset has fewer than 40 unique GO terms. Skipping.")
    return(invisible(NULL))
  }

  message("Create simMatrix")
  simMatrix <- tryCatch({
    rrvgo::calculateSimMatrix(bar_data$go_id, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
  }, error = function(e) {
    message("Failed to calculate simMatrix: ", e$message)
    return(NULL)
  })

  if (is.null(simMatrix) || !is.matrix(simMatrix) || length(dim(simMatrix)) != 2) {
    message("Invalid simMatrix structure. Skipping dataset: ", plot_name)
    return(invisible(NULL))
  }

  message("Reduce semantic similarity")
  scores <- stats::setNames(-log10(bar_data$p.adjust), bar_data$go_id)
  reducedTerms <- tryCatch({
    rrvgo::reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
  }, error = function(e) {
    message("reduceSimMatrix() failed: ", e$message)
    return(NULL)
  })

  if (is.null(reducedTerms)) {
    message("No reduced terms could be computed.")
    return(invisible(NULL))
  }

  reducedTerms <- reducedTerms[reducedTerms$size != 0, ]
  save(reducedTerms, simMatrix,
       file = file.path(outdir, paste0("reducedTerms_simMatrix_", direction, "_", plot_name, ".RData")))

  message("Merge metadata")
  joined_data <- dplyr::left_join(bar_data, reducedTerms, by = "term")
  joined_data <- joined_data[!is.na(joined_data$parentTerm), ]
  save(joined_data,
       file = file.path(outdir, paste0("geneListAndGOterms_", direction, "_", plot_name, ".RData")))

  message("Save treemap plot")
  grDevices::png(file = file.path(figdir, paste0("GOTerm_", direction, "_Treemap_", plot_name, ".png")))
  rrvgo::treemapPlot(reducedTerms)
  grDevices::dev.off()

  df <- dplyr::summarise(dplyr::group_by(reducedTerms, parentTerm), count = dplyr::n())
  df <- dplyr::arrange(df, dplyr::desc(count))

  bar_plot <- ggplot2::ggplot(utils::head(df, 50),
                              ggplot2::aes(x = stats::reorder(parentTerm, count), y = count)) +
    ggplot2::geom_bar(stat = "identity", fill = if (direction == "UP") "firebrick" else "dodgerblue") +
    ggplot2::labs(x = "Parent Terms", y = "Number of similar terms") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"),
                   panel.background = ggplot2::element_rect(fill = "white"))

  ggplot2::ggsave(filename = file.path(figdir, paste0("GOTerm_", direction, "_ReducedBar_", plot_name, ".png")),
                  plot = bar_plot, width = 30, height = 25, units = "cm", dpi = 300)
  invisible(NULL)
}


# --------- Reduce_GOTerm (function) / REDUCE GO TERMS----------------------
#' Batch Reduce and Visualize GO Terms Across Multiple Datasets
#'
#' Applies semantic similarity reduction to GO enrichment results from a batch
#' of \code{.RData} files (both UP and DOWN regulated genes). For each dataset,
#' it invokes \code{\link{processGO}} to reduce redundancy, merge metadata, and
#' generate visualizations.
#'
#' This wrapper function automates analysis across multiple files following the
#' naming pattern \code{UP_GO_<name>.RData} or \code{DOWN_GO_<name>.RData},
#' typically produced by GO enrichment pipelines.
#'
#' @param indir Character. Input directory containing GO enrichment
#'   \code{.RData} files. Expected filenames follow the pattern
#'   \code{UP_GO_<name>.RData} and \code{DOWN_GO_<name>.RData}.
#' @param outdir Character. Output directory to save reduced GO term
#'   \code{.RData} files.
#' @param figdir Character. Directory to save output plots (treemaps and bar plots).
#'
#' @details For each file, the function:
#' \itemize{
#'   \item Loads the GO enrichment result.
#'   \item Filters invalid, duplicated, or missing GO IDs.
#'   \item Skips datasets with fewer than 40 unique GO terms.
#'   \item Computes semantic similarity with \code{rrvgo::calculateSimMatrix()}
#'     using the "Rel" method.
#'   \item Reduces terms via \code{rrvgo::reduceSimMatrix()} with a similarity
#'     threshold of 0.7.
#'   \item Saves both the reduced term object and merged gene-GO mapping.
#'   \item Generates treemap and bar plots showing parent GO categories.
#' }
#' Internal errors in similarity computation are caught with \code{tryCatch()}
#' to allow uninterrupted batch processing across datasets.
#'
#' Requires the Bioconductor packages \pkg{rrvgo} and \pkg{org.Hs.eg.db}.
#'
#' @return Invisibly, \code{NULL}. It generates and saves \code{.RData} and
#'   \code{.png} files in the specified output directories.
#'
#' @seealso \code{\link{processGO}} for analyzing individual files.
#'
#' @examples
#' \dontrun{
#' # Batch process all GO enrichment results in the directory
#' Reduce_GOTerm(
#'   indir = "data/enrichment/",
#'   outdir = "results/GO_reduced/",
#'   figdir = "figures/GO_terms/"
#' )
#' }
#'
#' @export
Reduce_GOTerm <- function(indir, outdir, figdir) {
  files_up   <- list.files(indir, pattern = "^UP_GO_(.+)\\.RData$", full.names = TRUE)
  files_down <- list.files(indir, pattern = "^DOWN_GO_(.+)\\.RData$", full.names = TRUE)

  # Process both UP and DOWN datasets using processGO()
  for (f in files_up) {
    processGO(f = f, direction = "UP", outdir = outdir, figdir = figdir)
  }
  for (f in files_down) {
    processGO(f = f, direction = "DOWN", outdir = outdir, figdir = figdir)
  }
  invisible(NULL)
}

utils::globalVariables(c("parentTerm", "count"))
