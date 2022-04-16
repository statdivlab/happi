#' expit: inverse of logit
#'
#' @param input vector
#'
#' @export
expit <- function(input) {
  exp(input) / (1 + exp(input))
}

