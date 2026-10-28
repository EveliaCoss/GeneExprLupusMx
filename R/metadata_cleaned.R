#' Sample metadata: group, cell type and clinical information
#'
#' Metadata of the 99 RNA-seq samples: 33 volunteers (10 healthy controls and
#' 23 women with systemic lupus erythematosus, SLE), each contributing three
#' cell types: monocytes (\code{monocyte}), monocyte-derived dendritic cells
#' (\code{moDC}) and tolerogenic dendritic cells (\code{tolDC}). Clinical
#' variables are repeated in the three samples of each volunteer.
#'
#' @format A data frame with 99 rows (samples) and 14 columns. All columns are
#'   factors.
#' \describe{
#'   \item{sample_ID}{Sample identifier: volunteer ID plus cell type suffix
#'     (\code{_0} = monocyte, \code{_1} = moDC, \code{_3} = tolDC; e.g.,
#'     \code{"QR011_1"}). Also used as row names.}
#'   \item{Group}{\code{Ctrl} (healthy control) or \code{SLE} (patient).}
#'   \item{GC}{Glucocorticoid treatment: \code{No} or \code{Si}.}
#'   \item{Dose}{Glucocorticoid dose in mg (\code{0}, \code{2.5}, \code{3},
#'     \code{5} or \code{8}); \code{0} for controls and untreated patients.}
#'   \item{Dose_category}{Glucocorticoid treatment coded as \code{0} (no) or
#'     \code{1} (yes); equivalent to \code{GC}.}
#'   \item{Age}{Age in years.}
#'   \item{Cell_type}{Cell type: \code{monocyte}, \code{moDC} or \code{tolDC}.}
#'   \item{Group_age}{Age group: \code{less30} (younger than 30 years) or
#'     \code{more30} (30 years or older).}
#'   \item{Nephritis}{Lupus nephritis: \code{No} or \code{Si}; \code{NA} in
#'     controls.}
#'   \item{Other_tx}{Treatment other than glucocorticoids (e.g., antimalarials,
#'     immunosuppressants): \code{No} or \code{Si}; \code{NA} in controls.}
#'   \item{Sledai_score}{SLEDAI disease activity score (0 to 62); \code{NA} in
#'     controls.}
#'   \item{Sledai_group}{SLEDAI activity group: \code{Without_activity}
#'     (score 0), \code{Low} (1-4), \code{Moderate} (6), \code{High} (10-19)
#'     or \code{Very_high} (27 or more); \code{NA} in controls.}
#'   \item{SLEDAI_Md}{Yes/no indicator from the original clinical records
#'     (\code{No}, \code{Si} or \code{Sí}); \code{NA} in controls.}
#'   \item{SLEDAI_category}{Disease activity category: \code{Zero} (SLEDAI
#'     score 0, and all controls), \code{No_activity} (score 1-3) or
#'     \code{Activity} (score 4 or more).}
#' }
#'
#' @details
#' \itemize{
#'   \item Rows are in the same order as the columns of the count matrices of
#'     \code{\link{filtered_txi}}.
#'   \item \code{Age}, \code{Dose} and \code{Sledai_score} are stored as
#'     factors; convert them with \code{as.numeric(as.character(x))} before
#'     doing calculations.
#'   \item The volunteer identifier (e.g., \code{"QR011"}) is the prefix of
#'     \code{sample_ID}, \code{sub("_.*", "", sample_ID)}; it links this table
#'     with \code{\link{volunteer_info}}.
#'   \item Known inconsistencies of the original records: \code{SLEDAI_Md} has
#'     the level \code{Si} written in two ways (\code{Si} and \code{Sí}); in
#'     \code{Sledai_group}, one patient with score 0 is classified as \code{Low}
#'     and one with score 10 as \code{High}. Some ages and treatments differ
#'     from \code{\link{volunteer_info}}, which comes from the recruitment table.
#' }
#'
#' @source Clinical records of the study, curated for the analysis.
#'
#' @seealso \code{\link{filtered_txi}}, \code{\link{volunteer_info}}
#'
#' @examples
#' head(metadata_cleaned)
#' table(metadata_cleaned$Group, metadata_cleaned$Cell_type)
#'
#' # Age as a number
#' summary(as.numeric(as.character(metadata_cleaned$Age)))
"metadata_cleaned"
