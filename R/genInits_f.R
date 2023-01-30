#' Function for generating initial starts for estimated f
#'
#' @param nn length of outcome nn 
#' @param nstarts number of starts
#' @param seed numeric seed for random initializations
#' @param outcome vector of outcome, binary 0 & 1 
#'
#' @return matrix of initializations
#'
#' @examples
#' outcome <- rbinom(10, 1,.5)
#' genInits_f(nn = 10, nstarts = 1, seed = 88, outcome = outcome)
#' @export
genInits_f <- function(nn, 
                     nstarts = 1, 
                     seed, 
                     outcome) {
  
  init_start <- rep(mean(outcome), nn) # fixed initial start parameters
  inits <- init_start
  
  if (nstarts > 1) { 
    set.seed(seed)
    for (i in 2:nstarts) {
      inits <- rbind(inits, rep(runif(1,0,1),length(init_start)))
      # make sure to set.seed() before running multiple starts
    } 
  }
  
  return(inits)
}
