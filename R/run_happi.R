#' Main function for happi, p=q=1; this script contains the modularized version of happi with correct implementation of log likelihood
#'
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1
#' @param covariate n x p matrix; this is the matrix for the primary predictor/covariate of interest
#' @param quality_var length-n vector; this is the quality variable vector, currently p = 1  TODO(turn into n x q matrix)
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param nstarts number of starts; Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param change_threshold algorithm will terminate early if the likelihood changes by this percentage or less for 5 iterations in a row for both the alternative and the null
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1; default is 0. Either a single value or a vector of length n.  
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
#' @examples 
#' data(TM7_data)
#' x_matrix <- model.matrix(~tongue, data = TM7_data)
#' happi_results <- happi (outcome = TM7_data$`Cellulase/cellobiase CelA1`,
#' covariate=x_matrix, 
#' quality_var=TM7_data$mean_coverage,
#' max_iterations=1000, 
#' change_threshold=0.1,
#' epsilon=0, 
#' nstarts = 1, 
#' spline_df = 3)
#' @export
happi <- function(outcome,
                  covariate,
                  quality_var,
                  max_iterations = 1000,
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
  # TODO(PT) reduce redundancies i.e. calculation of ll 
  
  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data
  
  nn <- length(outcome)
  
  stopifnot(nn == nrow(covariate) | nn == length(quality_var))
  
  pp <- ncol(covariate)
  
  if (length(epsilon) == 1) {
    epsilon_vec <- rep(epsilon, nn)
  } else if (length(epsilon == nn)) {
    epsilon_vec <- epsilon
  } else {
    error("epsilon should be a single number or a length n vector.")
  }
  
  #  if (ncol(covariate) > 2) warning("Amy hasn't properly checked that multiple covariates result in sensible output")
  #  if(h0_param != 2) warning("Amy hasn't properly checked that testing a different parameter results in sensible output")
  
  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[my_order]
  outcome <- outcome[my_order]
  covariate <- covariate[my_order, ]
  epsilon_vec <- epsilon_vec[my_order]
  
  if (pp == 1) {
    covariate <- matrix(covariate, ncol = 1, nrow = nn)
    covariate_null <- matrix(1, ncol = 1, nrow = nn)
  } else {
    covariate_null <- covariate[, -h0_param]
  } #TODO(PT) take in formula 
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=nn)
  
  # generate beta initial starts for alternative model
  inits <- happi::genInits(num_covariate = pp, nstarts = nstarts, seed = seed)
  # generate beta initial starts for null model   
  inits_null <- happi::genInits(num_covariate = pp - 1, nstarts = nstarts, seed = seed)
  # generate f initial starts for both alternative and null models   
  inits_f <- happi::genInits_f(nn = nn, nstarts = nstarts, seed = seed, outcome = outcome)
  
  bestOut <- NULL
  bestOut_null <- NULL
  
  ### Start for loop through multiple starts 
  
  # save information from each start
  starts_df <- data.frame(starts = 1:nstarts,
                          alt_ll = NA,
                          null_ll = NA)
  
  for (i in 1:nstarts) {
    
    ####################################
    ## Create matrices to hold results #
    ####################################
    my_estimates <- tibble::tibble("iteration" = 0:max_iterations,
                                   "epsilon" = ifelse(length(epsilon) == 1, epsilon, "multiple"),
                                   "loglik" = NA,  
                                   "loglik_nopenalty" = NA)
    
    my_estimates_null <- tibble::tibble("iteration" = 0:max_iterations,
                                        "epsilon" = ifelse(length(epsilon) == 1, epsilon, "multiple"),
                                        "loglik" = NA, 
                                        "loglik_nopenalty" = NA)
    
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
    
    ##################
    ## Input starts ##
    ##################
    my_estimated_beta[1, ] <- inits[i,]
    my_estimated_beta_null[1, ] <- inits_null[i,]
    
    my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
    my_fitted_xbeta_null[1, ] <- c(covariate_null %*% my_estimated_beta_null[1, ])
    
    my_estimated_f[1, ] <- inits_f[i,]
    my_estimated_f_null[1, ] <- inits_f[i,]
    # f-tilde = logit(f)
    my_estimated_ftilde[1, ] <- happi::logit(my_estimated_f[1, ])
    my_estimated_ftilde_null[1, ] <- happi::logit(my_estimated_f_null[1, ])
    
    
    my_estimated_p[1, ] <- happi::calculate_p(xbeta = my_fitted_xbeta[1, ],
                                              ff = my_estimated_f[1, ], 
                                              epsilon = epsilon_vec, 
                                              outcome = outcome)
    my_estimated_p_null[1, ] <- happi::calculate_p(xbeta = my_fitted_xbeta_null[1, ],
                                                   ff = my_estimated_f_null[1, ], 
                                                   epsilon = epsilon_vec, 
                                                   outcome = outcome)
    
    my_estimates[1, "loglik"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                          ff = my_estimated_f[1, ],
                                                          firth = firth, 
                                                          outcome = outcome, 
                                                          epsilon = epsilon_vec, 
                                                          covariate = covariate)
    my_estimates[1, "loglik_nopenalty"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                                    ff = my_estimated_f[1, ],
                                                                    firth = F, 
                                                                    outcome = outcome, 
                                                                    epsilon = epsilon_vec, 
                                                                    covariate = covariate)
    my_estimates_null[1, "loglik"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                               ff = my_estimated_f_null[1, ],
                                                               firth = firth, 
                                                               outcome = outcome, 
                                                               epsilon = epsilon_vec, 
                                                               covariate = covariate_null)
    my_estimates_null[1, "loglik_nopenalty"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                                         ff = my_estimated_f_null[1, ],
                                                                         firth = F, 
                                                                         outcome = outcome, 
                                                                         epsilon = epsilon_vec, 
                                                                         covariate = covariate_null)
    ###############################################################################################################
    ## E-M algorithm for penalized maximum likelihood estimation of our parameter estimates for alternative model #
    ###############################################################################################################
    mlout <- tryCatch(happi::run_em(outcome = outcome,
                                    quality_var = quality_var,
                                    max_iterations = max_iterations,
                                    min_iterations = min_iterations,
                                    change_threshold = change_threshold,
                                    epsilon = epsilon_vec,
                                    method = method,
                                    firth = firth,
                                    spline_df = spline_df, 
                                    nn = nn, 
                                    em_covariate = covariate,
                                    em_estimated_p = my_estimated_p,
                                    em_fitted_xbeta = my_fitted_xbeta, 
                                    em_estimated_f = my_estimated_f, 
                                    em_estimates = my_estimates,
                                    em_estimated_beta  = my_estimated_beta,
                                    em_estimated_basis_weights =  my_estimated_basis_weights,
                                    em_estimated_ftilde = my_estimated_ftilde), 
                      error = function(e) {cat("WARNING alternative model E-M at initial start row index", paste(i),":", conditionMessage(e),"\n")}) 
    
    ########################################################################################################
    ## E-M algorithm for penalized maximum likelihood estimation of our parameter estimates for null model #
    ########################################################################################################
    mlout_null <- tryCatch(happi::run_em(outcome = outcome,
                                         quality_var = quality_var,
                                         max_iterations = max_iterations,
                                         min_iterations = min_iterations,
                                         change_threshold = change_threshold,
                                         epsilon = epsilon_vec,
                                         method = method,
                                         firth = firth,
                                         spline_df = spline_df, 
                                         nn = nn, 
                                         em_covariate = covariate_null,
                                         em_estimated_p = my_estimated_p_null,
                                         em_fitted_xbeta = my_fitted_xbeta_null, 
                                         em_estimated_f = my_estimated_f_null, 
                                         em_estimates = my_estimates_null,
                                         em_estimated_beta  = my_estimated_beta_null,
                                         em_estimated_basis_weights =  my_estimated_basis_weights_null,
                                         em_estimated_ftilde = my_estimated_ftilde_null), 
                           error = function(e) {cat("WARNING null model E-M at initial start row index", paste(i),":", conditionMessage(e),"\n")}) # END WHILE loop that stops when convergence is met
    ############################################################
    # Evaluate for best set of results for alternative model ###
    ############################################################
    ## if we get a converged result then... 
    ## check if bestOut is NULL. 
    ## If NULL then replace with the converged result 
    ## If NOT NUlL then replace only if LL is better than what is currently in bestOut
    
    tryCatch(if(!is.na(tail(mlout$loglik$loglik[!is.na(mlout$loglik$loglik)],1))){ # check that not NA for alternative model
      starts_df[i, "alt_ll"] <- tail(mlout$loglik$loglik[!is.na(mlout$loglik$loglik)],1)
      if (is.null(bestOut)) {
        bestOut <-  mlout
      } else if (!is.null(bestOut)){
        if(tail(mlout$loglik$loglik[!is.na(mlout$loglik$loglik)],1) > tail(bestOut$loglik$loglik[!is.na(bestOut$loglik$loglik)],1)){
          bestOut <- mlout 
        }
      }
    }, error = function(e) {cat("WARNING alternative model did not converge at initial start row index", paste(i),":", conditionMessage(e),"\n")})
    
    #####################################################
    # Evaluate for best set of results for null model ###
    #####################################################
    ## if we get a converged result then... 
    ## check if bestOut_null model is NULL. 
    ## If bestOut_null is NULL then replace with the converged result 
    ## If bestOut_null is NOT NUlL then replace only if LL is better than what is currently in bestOut_null
    
    tryCatch(if(!is.na(utils::tail(mlout_null$loglik$loglik[!is.na(mlout_null$loglik$loglik)],1))){ # check that not NA for null model 
      starts_df[i, "null_ll"] <- tail(mlout_null$loglik$loglik[!is.na(mlout_null$loglik$loglik)],1)
      if (is.null(bestOut_null)) {
        bestOut_null <-  mlout_null
      } else if (!is.null(bestOut_null)){
        if(utils::tail(mlout_null$loglik$loglik[!is.na(mlout_null$loglik$loglik)],1) > tail(bestOut_null$loglik$loglik[!is.na(bestOut_null$loglik$loglik)],1)){
          bestOut_null <- mlout_null 
        }
      }
    }, error = function(e) {cat("WARNING null model did not converge at initial start row index", paste(i),":", conditionMessage(e),"\n")})
    
  } # END -- nstarts
  
  if (is.null(bestOut)) stop("Full model could not be optimized! Try increasing the number of nstarts or simplifying your model.")
  if (is.null(bestOut_null)) stop("Null model could not be optimized! Try increasing the number of nstarts or simplifying your model.")
  
  #############################################
  ### If alternative LL is smaller than null:
  ## restart estimation of alternative model using the null model 
  #############################################
  if (utils::tail(bestOut$loglik$loglik_nopenalty[!is.na(bestOut$loglik$loglik_nopenalty)],1) < 
      tail(bestOut_null$loglik$loglik_nopenalty[!is.na(bestOut_null$loglik$loglik_nopenalty)],1)) {
    
    message("Likelihood greater under null; restarting...")
    my_estimates <- tibble::tibble("iteration" = 0:max_iterations,
                                   "epsilon" = ifelse(length(epsilon) == 1, epsilon, "multiple"),
                                   "loglik" = NA, 
                                   "loglik_nopenalty" = NA)
    my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
    my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    
    if (pp == 1){
      my_estimated_beta[1, 1] <- 0
      
    } else if (pp == 2) {
      my_estimated_beta[1, 1] <- utils::tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[1] ## start at converged null
      my_estimated_beta[1, h0_param] <- 0
      
    } else if (pp == 3){
      my_estimated_beta[1, 1] <- utils::tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[1] ## start at converged null
      my_estimated_beta[1, h0_param] <- 0
      my_estimated_beta[1, 3] <- tail(bestOut_null$beta[!is.na(bestOut_null$beta[,1]),],1)[2] ## start at converged null
      
    } ## TO DO (PT): need to take in pp > 3 if needed 
    
    stopifnot(h0_param == 2)
    my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
    my_estimated_f[1, ] <- utils::tail(bestOut_null$f[!is.na(bestOut_null$f[,1]),],1)
    my_estimated_ftilde[1, ] <- happi::logit(my_estimated_f[1, ])
    my_estimated_p[1, ] <- happi::calculate_p(xbeta = my_fitted_xbeta[1, ],
                                              ff = my_estimated_f[1, ], 
                                              epsilon = epsilon_vec, 
                                              outcome = outcome)
    my_estimates[1, "loglik"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                          ff = my_estimated_f[1, ], 
                                                          outcome = outcome, 
                                                          epsilon = epsilon_vec, 
                                                          covariate = covariate)
    my_estimates[1, "loglik_nopenalty"] <- happi::incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                                    ff = my_estimated_f[1, ], 
                                                                    outcome = outcome, 
                                                                    epsilon = epsilon_vec, 
                                                                    covariate = covariate,
                                                                    firth = F)
    ## restart at null model if likelihood greater under null than alternative
    bestOut <- happi::run_em(outcome = outcome,
                             quality_var = quality_var,
                             max_iterations = max_iterations,
                             min_iterations = min_iterations,
                             change_threshold = change_threshold,
                             epsilon = epsilon_vec,
                             method = method,
                             firth = firth,
                             nn = nn, 
                             spline_df = spline_df, 
                             em_covariate = covariate,
                             em_estimated_p = my_estimated_p,
                             em_fitted_xbeta = my_fitted_xbeta, 
                             em_estimated_f = my_estimated_f, 
                             em_estimates = my_estimates,
                             em_estimated_beta  = my_estimated_beta,
                             em_estimated_basis_weights =  my_estimated_basis_weights,
                             em_estimated_ftilde = my_estimated_ftilde)
    
    
    if (utils::tail(bestOut$loglik$loglik[!is.na(bestOut$loglik$loglik)],1) < utils::tail(bestOut_null$loglik$loglik[!is.na(bestOut_null$loglik$loglik)],1)) {
      message("Restarting to estimate beta_alt didn't work. Penalized likelihood is still greater under the null than alt. pvalue = 1")
      # message(paste("Had not converged after", tt_restart - 1, "iterations; LL % change:", round(pct_change_llks, 3)))
    }
    # if (tail(bestOut$loglik$loglik_nopenalty[!is.na(bestOut$loglik$loglik_nopenalty)],1) < tail(bestOut_null$loglik$loglik_nopenalty[!is.na(bestOut_null$loglik$loglik_nopenalty)],1)) {
    #    message("Weird! Nonpenalized Likelihood is also greater under the null than alt. pvalue = 1")
    #  } else {
    #    message("Phew! Nonpenalized Likelihood is not greater under the null than alt.")
    #  }
  } ## End restart if likelihood is greater under the null
  
  ###############################
  ### Output the best estimates## 
  ###############################
  my_estimates$loglik <- bestOut$loglik$loglik
  my_estimates$loglik_nopenalty <- bestOut$loglik$loglik_nopenalty
  my_estimates$iteration <- bestOut$loglik$iteration
  
  my_estimates$loglik_null <- bestOut_null$loglik$loglik
  my_estimates$loglik_null_nopenalty <- bestOut_null$loglik$loglik_nopenalty
  my_estimates$iteration_null <- bestOut_null$loglik$iteration
  
  my_estimated_beta <- bestOut$beta
  my_estimated_beta_null <- bestOut_null$beta
  my_estimated_f <- bestOut$f
  my_estimated_f_null <- bestOut_null$f
  my_estimated_basis_weights <- bestOut$basis_weights
  my_estimated_basis_weights_null <- bestOut_null$basis_weights
  my_estimated_p <- bestOut$p
  my_estimated_p_null <- bestOut_null$p
  
  ########################
  # LRT based on best LLs 
  ########################
  
  my_estimates$LRT <- 2*(utils::tail(my_estimates$loglik[!is.na(my_estimates$loglik)],1) - utils::tail(my_estimates$loglik_null[!is.na(my_estimates$loglik_null)],1))
  my_estimates$pvalue <- 1 - pchisq(my_estimates$LRT, df=1)
  
  my_estimates$LRT_nopenalty <- 2*(utils::tail(my_estimates$loglik_nopenalty[!is.na(my_estimates$loglik_nopenalty)],1) - utils::tail(my_estimates$loglik_null_nopenalty[!is.na(my_estimates$loglik_null_nopenalty)],1))
  my_estimates$pvalue_nopenalty <- 1 - pchisq(my_estimates$LRT_nopenalty, df=1)
  
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
              "covariate" = covariate,
              "starts_df" = starts_df))
  
  
} 
