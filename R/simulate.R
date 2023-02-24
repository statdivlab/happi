#' Internal simulate function for parametric bootstrap 
#'
#' @param happi_results_sims an object of happi results used for simulating new data
#' @param h0_param the column index in covariate that has beta=zero under the null
#'
#' @return a simulated dataframe using estimated null parameters
#' @export
#'
simulate_b <- function(happi_results_sims, 
                       h0_param = 2){
  my_data <- happi_results_sims
  
  beta_null <- as.vector(tail(my_data$beta_null[!is.na(my_data$beta_null[,1]),],1))
  f_null <- as.vector(tail(my_data$f_null[!is.na(my_data$f_null[,1]),],1))
  quality_var <- my_data$quality_var
  
  covariate <- my_data$covariate
  covariate_null <- covariate[, -h0_param]
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=length(covariate_null))
  
  lambda_props <- expit(covariate_null %*% beta_null)
  colnames(lambda_props) <- c("lambda_props")
  my_table <- cbind(f_null,quality_var,lambda_props,covariate) %>% as.data.frame()
  # Draw from a Bernoulli distribution with p = lambda_props 
  my_table$lambda <- sapply(my_table$lambda_props, function(x) rbinom(1, 1, x))
  # Probability of Y is 
  my_table$Y_prop <- sapply(my_table$f_null, function(x) expit(x))
  my_table$Y <- sapply(my_table$Y_prop, function(x) rbinom(1, 1, x))
  
  return(my_table)
}