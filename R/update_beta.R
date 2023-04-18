#' Function for updating betas 
#'
#' @param probs outcomes p; probability between 0 and 1 
#' @param covariate n x p matrix; this is the matrix for the primary predictor/covariate of interest
#' @param firth use firth penalty? Default is TRUE.
#'
#' @return vector of estimated betas
#'
#' @export
update_beta <- function(probs, covariate, firth = TRUE) {
  
  if (!firth) {
    # glm(probs ~ covariate - 1, family=binomial)$coef
    
    ## prevents warnings about `non-integer #successes in a binomial glm!`
    ## doesn't alter coefficient estimates compared to "binomial", only std errors, which we don't use
    coefs <- glm(probs ~ covariate - 1, family= quasibinomial)$coef
  } else {
    ### code that uses logistf (which seems not to have worked great...)
    #   if(all(covariate ==1)){ ## if covariate is just intercept, fit directly
    #     coefs <- logistf(probs~1)$coef
    #   } else{
    #   coefs <- logistf(probs ~ covariate - 1)$coef
    ### directly optimize penalized log likelihood:
    penalized_ll <-
      function(b){
        fitted_logits <- as.numeric(covariate%*%matrix(b,ncol = 1))
        pll <- sum(probs*fitted_logits - log(1 + exp(fitted_logits))) +
          0.5 * logdet(t(covariate) %*% diag(
            expit(fitted_logits) * (1 - expit(fitted_logits))) %*% covariate)
        return(-1*pll)
      }
    return(optim(rep(0,ncol(covariate)),penalized_ll, method = "L-BFGS-B")$par)
  }
}

