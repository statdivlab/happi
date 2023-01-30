#' Function for calculating incomplete log likelihood 
#'
#' @param xbeta my fitted betas matrix multiplied by covariate 
#' @param ff my estimated f's 
#' @param firth use Firth penalty? default: TRUE 
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1 
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param covariate n x p matrix; this is the matrix for the primary predictor/covariate of interest
#'
#' @return incomplete log-likelihood value
#'
#' @examples
#' data("TM7_data")
#' # create design matrix
#' covariate <- model.matrix(~tongue, data = TM7_data) 
#' outcome <- TM7_data$`Ribosomal protein L27`
#' beta <- c(0.1,0.03)
#' my_xbeta <- covariate %*% beta
#' nn <- length(outcome) 
#' my_ff <- rep(mean(outcome), nn)
#'
#' my_incomplete_LL <- incomplete_loglik(xbeta = my_xbeta, 
#' ff = my_ff, 
#' firth = TRUE, 
#' outcome = outcome, 
#' epsilon = 0, 
#' covariate = covariate)
#' 
#' @export
incomplete_loglik <- function(xbeta, 
                              ff, 
                              firth = TRUE, 
                              outcome = outcome, 
                              epsilon = epsilon, 
                              covariate = covariate) {
  
  prob_lambda <- expit(xbeta)
  
  if (!firth) {
    sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
               (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
  } else {
    ll <-  sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
                      (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
    penalty <- 0.5*msos::logdet(t(covariate) %*% diag(as.numeric(prob_lambda)*(
      1 - as.numeric(prob_lambda))) %*% covariate)
    return(ll + penalty)
  }
}
