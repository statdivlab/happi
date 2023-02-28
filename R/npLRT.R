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
#' @param seed numeric number to set seed for random multiple starts
#'
#' @return An object with npLRT pvalues 
#' @export
#'
npLRT <- function(happi_out,
                  max_iterations = 1000,
                  min_iterations = 15,
                  h0_param = 2,
                  change_threshold = 0.1,
                  epsilon = 0,
                  method = "splines",
                  firth = T,
                  spline_df = 4, 
                  nstarts = 1, 
                  seed = 8,
                  P = 1000){
  dif <- vector(length = P)
  set.seed(seed) # for reproducibility when we use the function sample() below
  dat <- happi_out
  random_rows <-sample(nrow(my_results$covariate))
  my_results$covariate[random_rows,]
  
  for (i in 1:length(dif)){
    dat$group <- sample(my_results$covariate) # shuffle the group labels
    # Compute the means for each group and then take the difference and store it
    dif[i] <- diff(tapply(X = dat$extra, 
                          INDEX = dat$group, 
                          FUN = mean))
  }
}