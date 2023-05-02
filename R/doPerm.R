#' Function to run a permutation iteration

#' @param happi_results_perm  a \code{happi} results object for permuting 
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param method method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by 
#' argument spline_df; "isotone" for isotonic regression fit
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit; default 3 
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param change_threshold algorithm will terminate early if the likelihood changes by this percentage or less for 
#' 5 iterations in a row for both the alternative and the null
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' 
#' @return test statistic for one permutation iteration 
#' @export
#'

doPerm <- function(happi_results_perm, 
                   h0_param = 2, 
                   method, 
                   firth, 
                   spline_df, 
                   max_iterations = 1000,
                   min_iterations = 15, 
                   change_threshold = 0.1, 
                   epsilon){
  
  my_results <- happi_results_perm
  
  random_rows <- sample(nrow(my_results$covariate))
  
  my_perm_covariate <- my_results$covariate[random_rows,] # shuffle the group labels
  
  happi_perm_out <- happi::happi(outcome=my_results$outcome, 
                    covariate=my_perm_covariate, 
                    quality_var=my_results$quality_var,
                    h0_param = h0_param, 
                    method = method, 
                    firth = firth, 
                    spline_df = spline_df, 
                    max_iterations = max_iterations,
                    min_iterations = min_iterations, 
                    change_threshold = change_threshold, 
                    epsilon = epsilon)
  
  LL_alt_nopenalty <- tail(happi_perm_out$loglik$loglik_nopenalty[!is.na(happi_perm_out$loglik$loglik_nopenalty)],1)
  LL_null_nopenalty <- tail(happi_perm_out$loglik$loglik_null_nopenalty[!is.na(happi_perm_out$loglik$loglik_null_nopenalty)],1)
  
  test.stat <- 2 * abs(LL_alt_nopenalty - LL_null_nopenalty)
  
  return(test.stat)
  
}