#' Test function running happi, p=q=1
#'
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1
#' @param covariate n x p matrix; this is the matrix for the primary predictor/covariate of interest
#' @param quality_var length-n vector; this is the quality variable vector, currently p = 1  TODO(turn into n x q matrix)
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param nstarts number of starts; Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param change_threshold algorithm will terminate early if the likelihood changes by this percentage or less for 5 iterations in a row for both the alternative and the null
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param method method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by 
#' argument spline_df; "isotone" for isotonic regression fit
#' @param random_starts whether to pick the starting values of beta's randomly. Defaults to FALSE.
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit; default 3 
#' @param seed numeric number to set seed for random multiple starts
#'
#' @importFrom msos logdet
#' @importFrom utils tail
#' @import tibble
#' @import stats
#' @import splines2
#' @importFrom isotone activeSet fSolver
#' @importFrom logistf logistf
#'
#' @return An object of class \code{happi}.
#'
#' @export
happi_test <- function(outcome,
                  covariate,
                  quality_var,
                  max_iterations = 50,
                  min_iterations = 15,
                  h0_param = 2,
                  change_threshold = 0.05,
                  epsilon = 0,
                  method = "splines",
                  random_starts = FALSE,
                  firth = TRUE,
                  spline_df = 3, 
                  nstarts = 1, 
                  seed = 13
) {
  
  # TODO(PT) take in formula
  
  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data
  
  nn <- length(outcome)
  
  stopifnot(nn == nrow(covariate) | nn == length(quality_var))
  
  pp <- ncol(covariate)
  
#  if (ncol(covariate) > 2) warning("Amy hasn't properly checked that multiple covariates result in sensible output")
#  if(h0_param != 2) warning("Amy hasn't properly checked that testing a different parameter results in sensible output")
  
  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[my_order]
  outcome <- outcome[my_order]
  covariate <- covariate[my_order, ]
  
  if (pp == 1) {
    covariate <- matrix(covariate, ncol = 1, nrow = nn)
    covariate_null <- matrix(1, ncol = 1, nrow = nn)
  } else {
    covariate_null <- covariate[, -h0_param]
  }
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=nn)
  
    inits <- genInits(num_covariate = pp, nstarts = nstarts, seed = seed)
    
    inits_null <- genInits(num_covariate = pp - 1, nstarts = nstarts, seed = seed)
    
    inits_f <- genInits_f(nn = nn, nstarts = nstarts, seed = seed, outcome = outcome)
    
    bestOut <- NULL
    bestOut_null <- NULL
    ### Optimization of the multiple starts in our penalized log likelihood 
    
    for (i in 1:nstarts) {
    my_estimates <- tibble("iteration" = 0:max_iterations,
                           "epsilon" = epsilon,
                           "loglik" = NA)
    my_estimates_null <- tibble("iteration" = 0:max_iterations,
                                "epsilon" = epsilon,
                                "loglik" = NA)
    
    my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
    my_estimated_beta_null <- matrix(NA, nrow = max_iterations + 1, ncol = max(1, pp - 1))
    my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_fitted_xbeta_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_f_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_ftilde_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_p_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_basis_weights <- matrix(NA, nrow = max_iterations + 1, ncol = spline_df + 1)
    my_estimated_basis_weights_null <- matrix(NA, nrow = max_iterations + 1, ncol = spline_df + 1)
    
    my_estimated_beta[1, ] <- inits[i,]
    my_estimated_beta_null[1, ] <- inits_null[i,]
    
    my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
    my_fitted_xbeta_null[1, ] <- c(covariate_null %*% my_estimated_beta_null[1, ])
    
    my_estimated_f[1, ] <- inits_f[i,]
    my_estimated_f_null[1, ] <- inits_f[i,]
    # f-tilde = logit(f)
    my_estimated_ftilde[1, ] <- logit(my_estimated_f[1, ])
    my_estimated_ftilde_null[1, ] <- logit(my_estimated_f_null[1, ])
    
    
    my_estimated_p[1, ] <- calculate_p(xbeta = my_fitted_xbeta[1, ],
                                       ff = my_estimated_f[1, ], 
                                       epsilon = epsilon, 
                                       outcome = outcome)
    my_estimated_p_null[1, ] <- calculate_p(xbeta = my_fitted_xbeta_null[1, ],
                                            ff = my_estimated_f_null[1, ], 
                                            epsilon = epsilon, 
                                            outcome = outcome)
    
    my_estimates[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                   ff = my_estimated_f[1, ],
                                                   firth = firth, 
                                                   outcome = outcome, 
                                                   epsilon = epsilon, 
                                                   covariate = covariate)
    
    my_estimates_null[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                        ff = my_estimated_f_null[1, ],
                                                        firth = firth, 
                                                        outcome = outcome, 
                                                        epsilon = epsilon, 
                                                        covariate = covariate_null)

    mlout <- tryCatch(run_em(outcome = outcome,
                        quality_var = quality_var,
                        max_iterations = max_iterations,
                        min_iterations = min_iterations,
                        change_threshold = change_threshold,
                        epsilon = epsilon,
                        method = method,
                        firth = firth,
                        spline_df = spline_df, 
                        em_covariate = covariate,
                        em_estimated_p = my_estimated_p,
                        em_fitted_xbeta = my_fitted_xbeta, 
                        em_estimated_f = my_estimated_f, 
                        em_estimates = my_estimates,
                        em_estimated_beta  = my_estimated_beta,
                        em_estimated_basis_weights =  my_estimated_basis_weights,
                        em_estimated_ftilde = my_estimated_ftilde), 
                        error = function(e) {cat("WARNING alternative model E-M at initial start row index", paste(i),":", conditionMessage(e),"\n")}) 
    
    mlout_null <- tryCatch(run_em(outcome = outcome,
                             quality_var = quality_var,
                             max_iterations = max_iterations,
                             min_iterations = min_iterations,
                             change_threshold = change_threshold,
                             epsilon = epsilon,
                             method = method,
                             firth = firth,
                             spline_df = spline_df, 
                             em_covariate = covariate_null,
                             em_estimated_p = my_estimated_p_null,
                             em_fitted_xbeta = my_fitted_xbeta_null, 
                             em_estimated_f = my_estimated_f_null, 
                             em_estimates = my_estimates_null,
                             em_estimated_beta  = my_estimated_beta_null,
                             em_estimated_basis_weights =  my_estimated_basis_weights_null,
                             em_estimated_ftilde = my_estimated_ftilde_null), 
                             error = function(e) {cat("WARNING null model E-M at initial start row index", paste(i),":", conditionMessage(e),"\n")}) # END WHILE loop that stops when convergence is met

    ## if we get a valid result then... check if bestOut is NULL. If NULL then replace with the valid result 
    
    tryCatch(if(!is.na(tail(mlout$loglik$loglik[!is.na(mlout$loglik$loglik)],1))){ # check that not NA
      if (is.null(bestOut)) {
        bestOut <-  mlout
      } else if (!is.null(bestOut)){
        if(tail(mlout$loglik$loglik[!is.na(mlout$loglik$loglik)],1) > tail(bestOut$loglik$loglik[!is.na(bestOut$loglik$loglik)],1)){
          bestOut <- mlout 
        }
      }
    },error = function(e) {cat("WARNING alternative model did not converge at initial start row index", paste(i),":", conditionMessage(e),"\n")})
    
    tryCatch(if(!is.na(tail(mlout_null$loglik$loglik[!is.na(mlout_null$loglik$loglik)],1))){ # check that not NA
      if (is.null(bestOut_null)) {
        bestOut_null <-  mlout_null
      } else if (!is.null(bestOut_null)){
        if(tail(mlout_null$loglik$loglik[!is.na(mlout_null$loglik$loglik)],1) > tail(bestOut_null$loglik$loglik[!is.na(bestOut_null$loglik$loglik)],1)){
          bestOut_null <- mlout_null 
        }
      }
    }, error = function(e) {cat("WARNING null model did not converge at initial start row index", paste(i),":", conditionMessage(e),"\n")})
    
  } # END -- nstarts
  
  if (is.null(bestOut)) stop("Full model could not be optimized! Try increasing the number of nstarts or simplifying your model.")
  if (is.null(bestOut_null)) stop("Null model could not be optimized! Try increasing the number of nstarts or simplifying your model.")
  
  
##### If alternative LL is smaller than null:       
     if (tail(bestOut$loglik$loglik[!is.na(bestOut$loglik$loglik)],1) < tail(bestOut_null$loglik$loglik[!is.na(bestOut_null$loglik$loglik)],1)) {
        
        message("Likelihood greater under null; restarting...")
        my_estimates <- tibble("iteration" = 0:max_iterations,
                              "epsilon" = epsilon,
                              "loglik" = NA)
        my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
        my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
        my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
        my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
        my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    
        if (pp == 1){
          my_estimated_beta[1, 1] <- 0
          
        } else if (pp == 2) {
          my_estimated_beta[1, 1] <- tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[1] ## start at converged null
          my_estimated_beta[1, h0_param] <- 0
          
        } else if (pp == 3){
          my_estimated_beta[1, 1] <- tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[1] ## start at converged null
          my_estimated_beta[1, h0_param] <- 0
          my_estimated_beta[1, 3] <- tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[2] ## start at converged null
          
        } ## TO DO PT 
        
        stopifnot(h0_param == 2)
        my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
        my_estimated_f[1, ] <- tail(bestOut_null$f[!is.na(bestOut_null$f[,1]),],1)
        my_estimated_ftilde[1, ] <- logit(my_estimated_f[1, ])
        my_estimated_p[1, ] <- calculate_p(xbeta = my_fitted_xbeta[1, ],
                                           ff = my_estimated_f[1, ], 
                                           epsilon = epsilon, 
                                           outcome = outcome)
        my_estimates[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                       ff = my_estimated_f[1, ], 
                                                       outcome = outcome, 
                                                       epsilon = epsilon, 
                                                       covariate = covariate)
        
        ## restart at null model if likelihood greater under null than alternative
        bestOut <- run_em(outcome = outcome,
                            quality_var = quality_var,
                            max_iterations = max_iterations,
                            min_iterations = min_iterations,
                            change_threshold = change_threshold,
                            epsilon = epsilon,
                            method = method,
                            firth = firth,
                            spline_df = spline_df, 
                            em_covariate = covariate,
                            em_estimated_p = my_estimated_p,
                            em_fitted_xbeta = my_fitted_xbeta, 
                            em_estimated_f = my_estimated_f, 
                            em_estimates = my_estimates,
                            em_estimated_beta  = my_estimated_beta,
                            em_estimated_basis_weights =  my_estimated_basis_weights,
                            em_estimated_ftilde = my_estimated_ftilde)
        

        if (tail(bestOut$loglik$loglik[!is.na(bestOut$loglik$loglik)],1) < tail(bestOut_null$loglik$loglik[!is.na(bestOut_null$loglik$loglik)],1)) {
          message("Restarting to estimate beta_alt didn't work. Likelihood is still greater under the null than alt.")
          # message(paste("Had not converged after", tt_restart - 1, "iterations; LL % change:", round(pct_change_llks, 3)))
        }
} ## End restart 
    ### Return the best estimates 
    
    my_estimates <- bestOut$loglik
    my_estimates$loglik_null <- bestOut_null$loglik$loglik
    my_estimates$iteration_null <- bestOut_null$loglik$iteration
    
    my_estimated_beta <- bestOut$beta
    my_estimated_beta_null <- bestOut_null$beta
    my_estimated_f <- bestOut$f
    my_estimated_f_null <- bestOut_null$f
    my_estimated_basis_weights <- bestOut$basis_weights
    my_estimated_basis_weights_null <- bestOut_null$basis_weights
    my_estimated_p <- bestOut$p
    my_estimated_p_null <- bestOut_null$p
    
    my_estimates$LRT <- 2*(tail(my_estimates$loglik[!is.na(my_estimates$loglik)],1) - tail(my_estimates$loglik_null[!is.na(my_estimates$loglik_null)],1))
    my_estimates$pvalue <- 1 - pchisq(my_estimates$LRT, df=1)

  return(list("loglik" = my_estimates,
              "beta" = my_estimated_beta,
              "beta_null" = my_estimated_beta_null,
              "f" = my_estimated_f,
              "f_null" = my_estimated_f_null,
              "basis_weights" =  my_estimated_basis_weights,
              "basis_weights_null" =  my_estimated_basis_weights_null,
              "p" = my_estimated_p,
              "p_null" = my_estimated_p_null,
              "quality_var" = quality_var,
              "outcome" = outcome,
              "covariate" = covariate))
  
  
} 

