#' parametric boostrap LRT implementation 
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
#' @return An object with pbLRT pvalues 
#' @export
#'
pbLRT <- function(happi_results_object, 
                  firth = TRUE, 
                  spline_df = 3, 
                  max_iterations = 1000, 
                  change_threshold = 0.01, 
                  epsilon = 0, 
                  nstarts = 1, 
                  h0_param = 2){
  
}