#' is.integer0 
#'
#' @description Function, that checks if certain integer is 0
#'
#' @return TRUE if 0, FALSE if not
#'
#' @noRd
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
#' %>%
#'
#' @description Pipe symbol from magrittr package
#'
#' @return %>%
#'
#' @noRd

