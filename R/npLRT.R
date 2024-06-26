#' Non parametric LRT aka permutation test
#'
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit; default 3 
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param change_threshold algorithm will terminate early if the likelihood changes
#' by this percentage or less for 5 iterations in a row for both the alternative and the null
#' @param epsilon probability of observing a gene when it should be absent; probability between 0 and 1
#' @param nstarts number of starts; Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param happi_out a happi results object 
#' @param P number of permutation iterations
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param method method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by 
#' argument spline_df; "isotone" for isotonic regression fit
#'
#' @return An object with npLRT pvalues 
#' @export
npLRT <- function(happi_out,
                  max_iterations = 1000,
                  min_iterations = 15,
                  h0_param = 2,
                  change_threshold,
                  epsilon,
                  method,
                  firth,
                  spline_df, 
                  nstarts = 1,
                  P = 1000){
  LL_alt <- utils::tail(happi_out$loglik$loglik_nopenalty[!is.na(happi_out$loglik$loglik_nopenalty)],1)
  LL_null<- utils::tail(happi_out$loglik$loglik_null_nopenalty[!is.na(happi_out$loglik$loglik_null_nopenalty)],1)
  
  # Compute LRT test statistic using data and estimated parameters of alternative and null models 
  t.observed <- 2 * (LL_alt-LL_null)    
  
  my_results_obj <- happi_out 
  PERM <- rep(NA, P)

  for (i in 1:P){
    # Compute the model 
    check_test_stat <- happi::doPerm(happi_results_perm = my_results_obj, 
                              spline_df = spline_df, 
                              max_iterations = max_iterations, 
                              h0_param = h0_param, 
                              epsilon = epsilon, 
                              min_iterations = min_iterations, 
                              change_threshold = change_threshold, 
                              method = method, 
                              firth = firth,
                              nstarts = nstarts)
    
    if (is.numeric(check_test_stat)){
      PERM[i] <- check_test_stat
    } else {
      PERM[i] <- NA
    }
  
  } # end PERM loop 
  
  if (sum(is.na(PERM)) == P) {
    warning("Nonparametric permutation test failed because the model failed to converge for all permutations.")
    return(NA)
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, PERM)
  message(paste0("Out of ", P, " permutations, ", sum(!is.na(PERM)), " converged and contributed to the permutation test."))
  return(p.val)
}



