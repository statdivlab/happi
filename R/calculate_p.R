#' Calculating p function 
#'
#' @param xbeta my fitted betas matrix multiplied by covariate 
#' @param ff my estimated f's 
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1
#'  
#'
#' @return a vector of probabilities  
#'
#' @examples
#' data("TM7_data")
#' outcome <- TM7_data$`Ribosomal protein L27`
#' # create design matrix
#' covariate <- model.matrix(~tongue, data = TM7_data) 
#' beta <- c(0.1,0.03)
#' my_xbeta <- covariate %*% beta
#' nn <- length(outcome)
#' my_ff <- rep(mean(outcome), nn)
#' epsilon <- 0 
#' 
#' my_calculated_p <- calculate_p(xbeta = my_xbeta, ff = my_ff, epsilon = 0, outcome = outcome)
#' @export
calculate_p <- function(xbeta, ff, epsilon, outcome) {
  pis <- ff * expit(xbeta) / (ff * expit(xbeta) + epsilon * (1 - expit(xbeta))) # for i: Y_i = 1
  pis_y0 <- (1 - ff) * expit(xbeta) / ((1 - ff) * expit(xbeta) + (1 - epsilon) * (1 - expit(xbeta)))
  pis[outcome == 0] <- pis_y0[outcome == 0]
  pis
}
