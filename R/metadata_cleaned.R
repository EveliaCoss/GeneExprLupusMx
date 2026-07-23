#' Cleaned metadata with patient and cell type annotations
#'
#' This dataset includes patient-level information along with cell type
#' annotations for monocytes (\code{mo}), monocyte-derived dendritic cells
#' (\code{moDC}), and tolerogenic dendritic cells (\code{tolDC}).
#'
#' @format A data frame with multiple rows (samples) and 13 columns:
#' \describe{
#'   \item{sample_ID}{Unique identifier for each sample.}
#'   \item{Group}{Disease group (e.g., SLE patient or control).}
#'   \item{GC}{Glucocorticoid treatment status.}
#'   \item{Dose}{Numeric dose of treatment.}
#'   \item{Dose_category}{Categorical dose classification.}
#'   \item{Age}{Age of the patient.}
#'   \item{Cell_type}{Annotated cell type (mo, moDC, tolDC).}
#'   \item{Group_age}{Combined group and age category.}
#'   \item{Nephritis}{Presence of nephritis (yes/no).}
#'   \item{Other_tx}{Other treatments received.}
#'   \item{Sledai_score}{SLEDAI disease activity score.}
#'   \item{Sledai_group}{Grouped SLEDAI classification.}
#'   \item{SLEDAI_Md}{Median SLEDAI score.}
#'   \item{SLEDAI_category}{Categorical SLEDAI classification.}
#' }
#'
#' @details This dataset is useful for linking transcriptomic data
#' with clinical and cell type annotations in systemic lupus erythematosus (SLE).
#'
#' @source Patient metadata curated and cleaned for downstream analysis.
#'
#' @examples
#' data(metadata_cleaned)
#' head(metadata_cleaned)
#' table(metadata_cleaned$Cell_type)
"metadata_cleaned"
