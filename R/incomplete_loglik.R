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
#' 
#' @export
incomplete_loglik <- function(xbeta, 
                              ff, 
                              firth = TRUE, 
                              outcome = outcome, 
                              epsilon = epsilon, 
                              covariate) {
  
  prob_lambda <- expit(xbeta)
  ## PT TODO: remove redundancy of the LL calculation 
  if (!firth) {
    sum(log( (1 - epsilon[outcome == 0]) * (1 - prob_lambda[outcome == 0]) +
               (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon[outcome == 1] * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
  } else {
    ll <-  sum(log( (1 - epsilon[outcome == 0]) * (1 - prob_lambda[outcome == 0]) +
                      (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon[outcome == 1] * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
    penalty <- 0.5*msos::logdet(t(covariate) %*% diag(as.numeric(prob_lambda)*(
      1 - as.numeric(prob_lambda))) %*% covariate)
    return(ll + penalty)
  }
}