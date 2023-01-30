#' Function for generating initial starts for estimated beta
#'
#' @param num_covariate number of covariates p (includes intercept)
#' @param nstarts number of starts
#' @param seed numeric seed for random initializations
#'
#' @return matrix of initializations
#'
#' @examples
#' genInits(num_covariate = 2, nstarts = 1, seed = 88)
#' @export
genInits <- function(num_covariate, 
                     nstarts = 1, 
                     seed) {
  
  init_start <- rep(0,num_covariate) # fixed initial start parameters
  inits <- rbind(init_start)

if (nstarts > 1) { 
    set.seed(seed)
    for (i in 2:nstarts) {
    inits <- rbind(inits, stats::rnorm(length(init_start), init_start, 15))
    # make sure to set.seed() before running multiple starts
    } 
  }
  
  return(inits)
}
