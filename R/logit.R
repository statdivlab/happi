#' logit for easy access
#'
#' @param input vector
#'
#' @export
logit <- function(input) {
  log(input / (1 - input))
}

