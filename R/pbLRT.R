#' parametric boostrap LRT implementation 
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
#' @param B number of bootstrap iterations
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param method method for estimating f. Defaults to "splines" which fits a monotone spline with df determined by 
#' argument spline_df; "isotone" for isotonic regression fit
#' @param seed numeric number to set seed for random multiple starts
#'
#' @return An object with pbLRT pvalues 
#' @export
#'
pbLRT <- function(happi_out,
                  max_iterations = 3000,
                  min_iterations = 15,
                  h0_param = 2,
                  change_threshold = 0.1,
                  epsilon = 0,
                  method = "splines",
                  firth = T,
                  spline_df = 4, 
                  nstarts = 1, 
                  seed = 88,
                  B = 1000){
# Step 1: Estimate parameters for alternative and null models 
# This function takes in penalized maximum likelihood estimates as input for 
# alternative and null models 
 LL_alt <- tail(happi_out$loglik$loglik_nopenalty[!is.na(happi_out$loglik$loglik_nopenalty)],1)
 LL_null<- tail(happi_out$loglik$loglik_null_nopenalty[!is.na(happi_out$loglik$loglik_null_nopenalty)],1)
 
# Compute LRT test statistic using data and estimated parameters of alternative and null models 
 t.observed <- 2 * (LL_alt-LL_null)    

 BOOT <- rep(NA, B)
 for (j in 1:B) {
   
   check_test_stat <- doBoot(happi_results_boot = happi_out)
   
   if (is.numeric(check_test_stat)){
     BOOT[j] <- check_test_stat
   } else {
     BOOT[j] <- NA
   }
   
 }
 perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
 p.val <- perc.rank(t.observed, BOOT)
 return(p.val)
}
 
 


