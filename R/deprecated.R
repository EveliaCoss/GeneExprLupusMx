#' Deprecated functions in GeneExprLupusMx
#'
#' These functions were renamed after version 0.1.0 to follow a single naming
#' style (\code{snake_case}). They still work, but show a warning and will be
#' removed in a future version.
#'
#' \itemize{
#'   \item \code{processGO()}: use \code{\link{process_go}()}.
#'   \item \code{Reduce_GOTerm()}: use \code{\link{reduce_go_terms}()}.
#' }
#'
#' @param ... Arguments passed to the new function.
#'
#' @return The value returned by the new function.
#'
#' @examples
#' # Use the new names:
#' # process_go()      instead of processGO()
#' # reduce_go_terms() instead of Reduce_GOTerm()
#' ?process_go
#'
#' @name GeneExprLupusMx-deprecated
#' @keywords internal
NULL

#' @rdname GeneExprLupusMx-deprecated
#' @export
processGO <- function(...) {
  .Deprecated("process_go", package = "GeneExprLupusMx")
  process_go(...)
}

#' @rdname GeneExprLupusMx-deprecated
#' @export
Reduce_GOTerm <- function(...) {
  .Deprecated("reduce_go_terms", package = "GeneExprLupusMx")
  reduce_go_terms(...)
}
