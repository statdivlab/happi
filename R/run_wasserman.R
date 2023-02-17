#' # Wasserman split test function 
#'
#' @param happi_results_object An object of class \code{happi}.
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit; default 3 
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param change_threshold algorithm will terminate early if the likelihood changes
#' by this percentage or less for 5 iterations in a row for both the alternative and the null
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param nstarts number of starts; Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param h0_param the column index in covariate that has beta=zero under the null
#'
#' @return A list of Wasserman pvalues 
#' @export
#'
run_wasserman <- function(happi_results_object, 
                          firth = TRUE, 
                          spline_df = 3, 
                          max_iterations = 1000, 
                          change_threshold = 0.01, 
                          epsilon = 0, 
                          nstarts = 1, 
                          h0_param = 2) {
  #take in the data 
  picked <- sample(seq_len(length(happi_results_object$outcome)), size = length(happi_results_object$outcome)/2, replace = FALSE)
  # split the data randomly into sets D0 and D1
  D1_results <- happi(outcome=happi_results_object$outcome[picked], 
                      covariate=happi_results_object$covariate[picked,], 
                      quality_var=happi_results_object$quality_var[picked],
                      method = "splines", 
                      firth = firth, 
                      spline_df = spline_df,
                      max_iterations=max_iterations, 
                      change_threshold= change_threshold, 
                      epsilon=epsilon, 
                      nstarts = nstarts)
  
  D0_results <- happi(outcome=happi_results_object$outcome[-picked], 
                      covariate=happi_results_object$covariate[-picked,], 
                      quality_var=happi_results_object$quality_var[-picked],
                      method="splines", 
                      firth=firth, 
                      spline_df= spline_df,
                      max_iterations=max_iterations, 
                      change_threshold=change_threshold, 
                      epsilon=epsilon, 
                      nstarts = nstarts)
  
  D1_outcome <- D1_results$outcome
  D1_quality_var <- D1_results$quality_var
  D1_covariate <- D1_results$covariate
  D1_covariate_null <- D1_covariate[, -h0_param]
  if (!is.matrix(D1_covariate_null)) D1_covariate_null <- matrix(D1_covariate_null, nrow=length(D1_covariate_null))
  
  D0_outcome <- D0_results$outcome
  D0_quality_var <- D0_results$quality_var
  D0_covariate <- D0_results$covariate
  D0_covariate_null <- D0_covariate[, -h0_param]
  if (!is.matrix(D0_covariate_null)) D0_covariate_null <- matrix(D0_covariate_null, nrow=length(D0_covariate_null))
  
  #################
  # Calculate U_n #
  #################
  D1_fitted_xbeta_D0_covariate<- c(D0_covariate %*% as.vector(tail(D1_results$beta[!is.na(D1_results$beta[,1]),],1)))
  D1_estimated_f <- as.vector(tail(D1_results$f[!is.na(D1_results$f[,1]),],1))
 
  LL_D0_theta_D1 <- incomplete_loglik(xbeta = D1_fitted_xbeta_D0_covariate,
                                       ff = D1_estimated_f,
                                       firth = T, 
                                       outcome = D0_outcome, 
                                       epsilon = 0, 
                                       covariate = D0_covariate)
  LL_D0_theta_D1_nopenalty <- incomplete_loglik(xbeta = D1_fitted_xbeta_D0_covariate,
                                      ff = D1_estimated_f,
                                      firth = F, 
                                      outcome = D0_outcome, 
                                      epsilon = 0, 
                                      covariate = D0_covariate)
  
  
  D1_fitted_xbeta_null_D0_covariate_null <- c(D0_covariate_null %*% as.vector(tail(D1_results$beta_null[!is.na(D1_results$beta_null[,1]),],1)))
  D1_estimated_f_null <- as.vector(tail(D1_results$f_null[!is.na(D1_results$f_null[,1]),],1))
  
  LL_D0_theta_D1_null <- incomplete_loglik(xbeta = D1_fitted_xbeta_null_D0_covariate_null,
                                     ff = D1_estimated_f_null, 
                                     firth = T,
                                     outcome = D0_outcome,
                                     epsilon = 0, 
                                     covariate = D0_covariate_null)
  LL_D0_theta_D1_null_nopenalty <- incomplete_loglik(xbeta = D1_fitted_xbeta_null_D0_covariate_null,
                                           ff = D1_estimated_f_null, 
                                           firth = F,
                                           outcome = D0_outcome,
                                           epsilon = 0, 
                                           covariate = D0_covariate_null)
  U_n <- exp(LL_D0_theta_D1-LL_D0_theta_D1_null)
  U_n_nopenalty <- exp(LL_D0_theta_D1_nopenalty-LL_D0_theta_D1_null_nopenalty)
  
  ######################
  # Calculate U_n_swap #
  ######################
  D0_fitted_xbeta_D1_covariate<- c(D1_covariate %*% as.vector(tail(D0_results$beta[!is.na(D0_results$beta[,1]),],1)))
  D0_estimated_f <- as.vector(tail(D0_results$f[!is.na(D0_results$f[,1]),],1))
  
  LL_D1_theta_D0 <- incomplete_loglik(xbeta = D0_fitted_xbeta_D1_covariate,
                                      ff = D0_estimated_f,
                                      firth = T, 
                                      outcome = D1_outcome, 
                                      epsilon = 0, 
                                      covariate = D1_covariate)
  LL_D1_theta_D0_nopenalty <- incomplete_loglik(xbeta = D0_fitted_xbeta_D1_covariate,
                                      ff = D0_estimated_f,
                                      firth = F, 
                                      outcome = D1_outcome, 
                                      epsilon = 0, 
                                      covariate = D1_covariate)
  
  D0_fitted_xbeta_null_D1_covariate_null <- c(D1_covariate_null %*% as.vector(tail(D0_results$beta_null[!is.na(D0_results$beta_null[,1]),],1)))
  D0_estimated_f_null <- as.vector(tail(D0_results$f_null[!is.na(D0_results$f_null[,1]),],1))
  
  LL_D1_theta_D0_null <- incomplete_loglik(xbeta = D0_fitted_xbeta_null_D1_covariate_null,
                                           ff = D0_estimated_f_null, 
                                           firth = T,
                                           outcome = D1_outcome,
                                           epsilon = 0, 
                                           covariate = D1_covariate_null)
  LL_D1_theta_D0_null_nopenalty <- incomplete_loglik(xbeta = D0_fitted_xbeta_null_D1_covariate_null,
                                           ff = D0_estimated_f_null, 
                                           firth = F,
                                           outcome = D1_outcome,
                                           epsilon = 0, 
                                           covariate = D1_covariate_null)
  U_n_swap <- exp(LL_D1_theta_D0-LL_D1_theta_D0_null)
  U_n_swap_nopenalty <- exp(LL_D1_theta_D0_nopenalty-LL_D1_theta_D0_null_nopenalty)
  
  pvalue <- 1/((U_n + U_n_swap)/2)
  pvalue_nopenalty <- 1/((U_n_nopenalty + U_n_swap_nopenalty)/2)

  return(list("wasserman_pvalue" = pvalue,
              "wassernan_pvalue_nopenalty" = pvalue_nopenalty))
}