#' Clinical information per volunteer
#'
#' Clinical and recruitment information with one row per volunteer
#' (10 healthy controls and 23 women with systemic lupus erythematosus).
#' Each volunteer contributed three samples (monocyte, moDC and tolDC) to
#' \code{metadata_cleaned}; both tables are linked by \code{Volunteer_ID}, the
#' prefix of \code{metadata_cleaned$sample_ID} (e.g. \code{"QR011_0"} -> \code{"QR011"}).
#'
#' @format A data frame with 33 rows and 16 columns:
#' \describe{
#'   \item{Volunteer_ID}{Volunteer identifier (e.g. \code{"QR011"}).}
#'   \item{Volunteer}{Volunteer number, in order of recruitment.}
#'   \item{Group}{\code{Ctrl} (healthy control) or \code{SLE} (patient).}
#'   \item{Age_years}{Age in years.}
#'   \item{Years_with_diagnosis}{Years since SLE diagnosis (\code{NA} in controls).}
#'   \item{Age_at_diagnosis_years}{Age at SLE diagnosis (\code{NA} in controls).}
#'   \item{SLEDAI}{SLEDAI disease activity: \code{1} = active (SLEDAI >= 6 points),
#'     \code{0} = no or low activity (0-5 points); \code{NA} in controls.}
#'   \item{Nephritis}{Lupus nephritis: \code{No} or \code{Yes}.}
#'   \item{Treatment_GCC}{Receives glucocorticoids (\code{TRUE}/\code{FALSE}).}
#'   \item{GCC_drug}{Glucocorticoid as recorded: \code{PDN} or \code{PDS}.}
#'   \item{GCC_dose_mg}{Glucocorticoid dose in mg (\code{NA} if not recorded).}
#'   \item{Treatment_AM}{Receives antimalarials.}
#'   \item{Treatment_MMF}{Receives mycophenolate mofetil.}
#'   \item{Treatment_AZT}{Receives azathioprine.}
#'   \item{Treatment_MTT}{Receives methotrexate.}
#'   \item{Time_of_recruitment}{Month of recruitment and sample collection
#'     (\code{Date}, first day of the month).}
#' }
#'
#' @details Treatment columns are \code{NA} for controls (not applicable).
#'   The cleaning script is in \code{data-raw/volunteer_info.R}.
#'
#' @source Hernández-Ledesma AL, Coss-Navarrete EL, et al. (2026), Front.
#'   Immunol. 17:1771959 (\doi{10.3389/fimmu.2026.1771959}); data and analysis
#'   code deposited in Zenodo (\doi{10.5281/zenodo.17742206}).
#'
#' @examples
#' data(volunteer_info)
#' head(volunteer_info)
#' table(volunteer_info$Group, format(volunteer_info$Time_of_recruitment, "%Y-%m"))
"volunteer_info"
