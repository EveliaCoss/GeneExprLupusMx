#' Gene-level RNA-seq quantification (tximport)
#'
#' Gene-level RNA-seq quantification of the 99 samples of the study
#' (33 volunteers x 3 cell types: monocytes, moDCs and tolDCs). Transcripts were
#' quantified with \pkg{Salmon} and summarized to genes with \pkg{tximport}.
#' Genes are identified by their gene symbol; duplicated entries in gene
#' biotypes were removed before creating this object.
#'
#' @format A list with the structure returned by \code{tximport::tximport()}:
#' \describe{
#'   \item{abundance}{Numeric matrix (21,718 genes x 99 samples) of gene
#'     abundances in transcripts per million (TPM).}
#'   \item{counts}{Numeric matrix (21,718 genes x 99 samples) of estimated
#'     read counts. Values are not integers because Salmon estimates them.}
#'   \item{length}{Numeric matrix (21,718 genes x 99 samples) with the average
#'     effective transcript length of each gene in each sample; \pkg{DESeq2}
#'     uses it as a normalization offset.}
#'   \item{countsFromAbundance}{Character, \code{"no"}: counts were not
#'     derived from abundances.}
#' }
#'
#' @details Row names are gene symbols (e.g., \code{"A1BG"}) and column names
#'   are sample identifiers (e.g., \code{"QR011_1"}), in the same order as the
#'   rows of \code{\link{metadata_cleaned}}. Use it with
#'   \code{DESeq2::DESeqDataSetFromTximport()}.
#'
#' @source Hernández-Ledesma AL, Coss-Navarrete EL, et al. (2026), Front.
#'   Immunol. 17:1771959 (\doi{10.3389/fimmu.2026.1771959}); analysis deposited
#'   in Zenodo (\doi{10.5281/zenodo.17742206}) and code in
#'   \url{https://github.com/NeuroGenomicsMX/SLE_Transcriptomic_analysis_tolerogenic}.
#'
#' @seealso \code{\link{metadata_cleaned}}, \code{\link{biotype_summary}}
#'
#' @examples
#' names(filtered_txi)
#' dim(filtered_txi$counts)
#' filtered_txi$counts[1:4, 1:4]
"filtered_txi"
