#' General E-M Algorithm 
#'
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1
#' @param em_covariate n x p matrix; this is the matrix for the primary predictor/covariate of interest
#' @param quality_var length-n vector; this is the quality variable vector, currently p = 1  TODO(turn into n x q matrix)
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param em_estimated_p estimated probablities
#' @param em_fitted_xbeta fitted betas
#' @param em_estimated_f estimated f's 
#' @param change_threshold algorithm will terminate early if the likelihood changes by this percentage or less for 5 iterations in a row for both th
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param method method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by 
#' argument spline_df; "isotone" for isotonic regression fit
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit; default 3 
#' @param em_estimates log likelihood estimates 
#' @param em_estimated_beta estimated betas
#' @param em_estimated_basis_weights estimated basis weights
#' @param em_estimated_ftilde estimated f_tilde aka logit(estimated_f)
#' @param nn length(outcome)
#'
#' @importFrom utils tail
#' 
#' @return An object of class \code{happi}.
#' @export
#'
 
run_em <- function(outcome = outcome,
                   quality_var = quality_var, 
                   change_threshold = change_threshold,
                   max_iterations = max_iterations,
                   min_iterations = min_iterations,
                   epsilon = epsilon,
                   method = method,
                   firth = firth,
                   spline_df = spline_df, 
                   nn = nn, 
                   em_covariate = NULL,
                   em_estimates = NULL,
                   em_estimated_beta  = NULL,
                   em_estimated_basis_weights =  NULL,
                   em_estimated_ftilde = NULL, 
                   em_estimated_p = NULL,
                   em_fitted_xbeta = NULL, 
                   em_estimated_f = NULL){

  tt <- 1
  keep_going <- TRUE
  
  while (tt <= max_iterations & keep_going) {
    
    tt <- tt + 1
    
    em_estimated_beta[tt, ] <- update_beta(probs=em_estimated_p[tt - 1, ],
                                           covariate=em_covariate,
                                           firth = firth)
    em_fitted_xbeta[tt, ] <- c(em_covariate %*% em_estimated_beta[tt, ])
    
    updated_f <- update_f(probs=em_estimated_p[tt - 1, ],
                          method = method,
                          spline_df = spline_df, 
                          outcome = outcome, 
                          nn = nn, 
                          quality_var = quality_var, 
                          tt = tt)
    em_estimated_ftilde[tt, ] <- updated_f$fitted_f_tilde
    em_estimated_basis_weights[tt, ] <- updated_f$basis_weights
    
    em_estimated_f[tt, ] <- expit(em_estimated_ftilde[tt, ])
    
    em_estimated_p[tt, ] <- calculate_p(xbeta=em_fitted_xbeta[tt, ], 
                                        ff=em_estimated_f[tt, ], 
                                        epsilon = epsilon, 
                                        outcome = outcome)
    
    em_estimates[tt, "loglik"] <- incomplete_loglik(xbeta = em_fitted_xbeta[tt, ],
                                                    ff = em_estimated_f[tt, ],
                                                    firth = firth, 
                                                    outcome = outcome, 
                                                    epsilon = epsilon, 
                                                    covariate = em_covariate)
    
    em_estimates[tt, "loglik_nopenalty"] <- incomplete_loglik(xbeta = em_fitted_xbeta[tt, ],
                                                    ff = em_estimated_f[tt, ],
                                                    firth = F, 
                                                    outcome = outcome, 
                                                    epsilon = epsilon, 
                                                    covariate = em_covariate)
    
    if ((tt > min_iterations) & (!is.na(em_estimates[tt, "loglik"]))) {
      
      pct_change_llks <- 100*abs((em_estimates[(tt - 4):tt, "loglik"] - em_estimates[(tt - 5):(tt - 1), "loglik"])/em_estimates[(tt - 5):(tt - 1), "loglik"])
      keep_going <- pct_change_llks > change_threshold # when this is all TRUE 
      
      if(all(!keep_going)){
        message(paste("Model converged after", tt, "iterations; LL % change:", round(tail(pct_change_llks,1), 3)))
        keep_going <- FALSE 
        check_updated_f <- warningcheck_update_f(probs=em_estimated_p[tt - 1, ],
                                                 method = method,
                                                 spline_df = spline_df, 
                                                 quality_var = quality_var, 
                                                 nn = nn, 
                                                 outcome = outcome, 
                                                 tt = tt)
        } else if((tt == max_iterations) & keep_going) {
          stop(paste("Model did not converge after", tt, "maximum number of iterations"))
         
        } else {
          keep_going <- TRUE 
        }
      
    } # END if - convergence of LL's 
    
  }
  
  return(list("loglik" = em_estimates,
              "beta" = em_estimated_beta,
              "f" = em_estimated_f,
              "basis_weights" =  em_estimated_basis_weights,
              "p" = em_estimated_p,
              "quality_var" = quality_var,
              "outcome" = outcome,
              "covariate" = em_covariate))
  
}
