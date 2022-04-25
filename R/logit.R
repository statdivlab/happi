#' logit: function for logit calculation
#'
#' @param input vector of value(s)  consisting of proportions `p` between 0 and 1
#'
#' @return Returns value(s) between negative infinity and infinity 
#' @export
#' 
#' @examples 
#' # taking the logit of any proportion between 0 and 1
#' logit(0.4)
#' # taking the logit of a vector of proportions 
#' logit(c(0.1, 0.2, 0.9))
logit <- function(input) {
  log(input / (1 - input))
}

