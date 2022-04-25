#' expit: inverse of logit
#'
#' @param input input a vector of value(s). These values can be between negative infinity to infinity.
#' 
#' @return Returns value(s) that are proportion(s) between 0 and 1 
#' @export
#' 
#' @examples 
#' # taking the expit of any value(s) between negative infinity and infinity & return a proportion 
#' expit(3)
#' # taking the expit of a vector of values 
#' expit(c(-0.1, 2, 1.5, -1))
expit <- function(input) {
  exp(input) / (1 + exp(input))
}

