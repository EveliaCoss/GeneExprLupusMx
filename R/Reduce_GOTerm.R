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
#' \donttest{
#' if (all(vapply(c("rrvgo", "org.Hs.eg.db", "dplyr", "ggplot2"),
#'                requireNamespace, logical(1), quietly = TRUE))) {
#'   # Enrichment table with the columns go_id, term and p.adjust
#'   go_bp <- go_celltype_results[go_celltype_results$query == "SLE_UP_moDC" &
#'                                  go_celltype_results$source == "GO:BP", ]
#'   go_bp <- head(go_bp[order(go_bp$p_value), ], 60)
#'   go_table <- data.frame(go_id = go_bp$term_id, term = go_bp$term_name,
#'                          p.adjust = go_bp$p_value)
#'
#'   go_dir <- file.path(tempdir(), "processGO_example")
#'   dir.create(go_dir, showWarnings = FALSE)
#'   f <- file.path(go_dir, "UP_GO_SLE_moDC.RData")
#'   save(go_table, file = f)
#'
#'   processGO(f, direction = "UP", outdir = go_dir, figdir = go_dir)
#'   list.files(go_dir)
#' }
#' }
#'
#' @export
processGO <- function(f, direction = c("UP", "DOWN"), outdir, figdir) {
  pkgs <- c("rrvgo", "org.Hs.eg.db", "dplyr", "ggplot2")
  if (!all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Install rrvgo, org.Hs.eg.db, dplyr and ggplot2 to use this function: ",
         "BiocManager::install(c('rrvgo', 'org.Hs.eg.db', 'dplyr', 'ggplot2'))")
  }
  direction <- match.arg(direction)
  plot_name <- gsub(paste0("^", direction, "_GO_(.+)\\.RData$"), "\\1", basename(f))

  message("\n------------------------------------")
  message("Dataset: ", plot_name)
  message("------------------------------------")

  bar_data <- read_go_file(f)
  if (is.null(bar_data)) {
    return(invisible(NULL))
  }
  reduced <- reduce_go_similarity(bar_data, plot_name)
  if (is.null(reduced)) {
    return(invisible(NULL))
  }
  save_reduced_go(bar_data, reduced, direction, plot_name, outdir)
  plot_reduced_go(reduced$reducedTerms, direction, plot_name, figdir)
  invisible(NULL)
}

# Carga un .RData con una sola tabla de enriquecimiento y quita los go_id
# vacíos o repetidos. Regresa NULL (con un mensaje) si no se puede analizar.
#' @noRd
read_go_file <- function(f) {
  env <- new.env() # entorno aparte, sin contaminar el resto
  load(f, envir = env)
  obj_names <- ls(env)
  if (length(obj_names) != 1) {
    stop("El archivo debe contener exactamente un objeto. Se encontraron: ",
         toString(obj_names))
  }

  bar_data <- env[[obj_names]]
  if (!("go_id" %in% colnames(bar_data))) {
    message("No 'go_id' column found. Skipping file: ", f)
    return(NULL)
  }

  bar_data <- bar_data[!is.na(bar_data$go_id) & !duplicated(bar_data$go_id), ]
  if (nrow(bar_data) < 40) {
    message("Dataset has fewer than 40 unique GO terms. Skipping.")
    return(NULL)
  }
  bar_data
}

# Agrupa los términos GO por similitud semántica (rrvgo). Regresa
# list(simMatrix, reducedTerms) o NULL (con un mensaje) si falla.
#' @noRd
reduce_go_similarity <- function(bar_data, plot_name) {
  message("Create simMatrix")
  simMatrix <- tryCatch(
    rrvgo::calculateSimMatrix(bar_data$go_id, orgdb = "org.Hs.eg.db",
                              ont = "BP", method = "Rel"),
    error = function(e) {
      message("Failed to calculate simMatrix: ", e$message)
      NULL
    }
  )
  if (!is.matrix(simMatrix)) {
    message("Invalid simMatrix structure. Skipping dataset: ", plot_name)
    return(NULL)
  }

  message("Reduce semantic similarity")
  scores <- stats::setNames(-log10(bar_data$p.adjust), bar_data$go_id)
  reducedTerms <- tryCatch(
    rrvgo::reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db"),
    error = function(e) {
      message("reduceSimMatrix() failed: ", e$message)
      NULL
    }
  )
  if (is.null(reducedTerms)) {
    message("No reduced terms could be computed.")
    return(NULL)
  }
  list(simMatrix = simMatrix, reducedTerms = reducedTerms[reducedTerms$size != 0, ])
}

# Guarda los términos reducidos y su unión con la tabla original
#' @noRd
save_reduced_go <- function(bar_data, reduced, direction, plot_name, outdir) {
  reducedTerms <- reduced$reducedTerms
  simMatrix <- reduced$simMatrix
  suffix <- paste0(direction, "_", plot_name, ".RData")
  save(reducedTerms, simMatrix,
       file = file.path(outdir, paste0("reducedTerms_simMatrix_", suffix)))

  message("Merge metadata")
  joined_data <- dplyr::left_join(bar_data, reducedTerms, by = "term")
  joined_data <- joined_data[!is.na(joined_data$parentTerm), ]
  save(joined_data, file = file.path(outdir, paste0("geneListAndGOterms_", suffix)))
}

# Guarda el treemap y la gráfica de barras de los términos padre
#' @noRd
plot_reduced_go <- function(reducedTerms, direction, plot_name, figdir) {
  message("Save treemap plot")
  save_treemap(reducedTerms,
               file.path(figdir, paste0("GOTerm_", direction, "_Treemap_", plot_name, ".png")))

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
}

# Dibuja el treemap de rrvgo en un PNG; el dispositivo se cierra aunque falle
#' @noRd
save_treemap <- function(reducedTerms, path) {
  grDevices::png(filename = path)
  on.exit(grDevices::dev.off(), add = TRUE)
  rrvgo::treemapPlot(reducedTerms)
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
#' \donttest{
#' if (all(vapply(c("rrvgo", "org.Hs.eg.db", "dplyr", "ggplot2"),
#'                requireNamespace, logical(1), quietly = TRUE))) {
#'   indir <- file.path(tempdir(), "Reduce_GOTerm_example")
#'   dir.create(indir, showWarnings = FALSE)
#'
#'   # One file per enrichment result, named UP_GO_<name>.RData or DOWN_GO_<name>.RData
#'   for (q in c("SLE_UP_moDC", "SLE_DOWN_moDC")) {
#'     go_bp <- go_celltype_results[go_celltype_results$query == q &
#'                                    go_celltype_results$source == "GO:BP", ]
#'     go_bp <- head(go_bp[order(go_bp$p_value), ], 60)
#'     go_table <- data.frame(go_id = go_bp$term_id, term = go_bp$term_name,
#'                            p.adjust = go_bp$p_value)
#'     direction <- if (grepl("_UP_", q, fixed = TRUE)) "UP" else "DOWN"
#'     save(go_table, file = file.path(indir, paste0(direction, "_GO_moDC.RData")))
#'   }
#'
#'   Reduce_GOTerm(indir = indir, outdir = tempdir(), figdir = tempdir())
#' }
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
