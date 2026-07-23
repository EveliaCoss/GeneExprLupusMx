#' Filtered transcript count data
#'
#' This dataset contains transcript-level count data obtained with \pkg{Salmon},
#' imported using \pkg{tximport}, and annotated with Ensembl gene identifiers.
#' Duplicate entries in gene biotypes were previously removed to ensure
#' consistency in downstream analyses.
#'
#' @format A list object with the following elements:
#' \describe{
#'   \item{counts}{Matrix of transcript counts.}
#'   \item{abundance}{Matrix of normalized transcript abundances.}
#'   \item{length}{Effective transcript lengths.}
#'   \item{metadata}{Data frame containing sample-level metadata.}
#' }
#'
#' @details This object is typically used for differential gene expression
#' analysis with \pkg{DESeq2}, after importing transcript quantifications
#' and mapping them to Ensembl identifiers.
#'
#' @source Processed internally from transcript quantifications generated with Salmon.
#'
#' @examples
#' data(filtered_txi)
#' names(filtered_txi)
#' head(filtered_txi$counts)
"filtered_txi"
