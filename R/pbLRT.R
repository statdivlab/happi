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
#'
#' @return An object with pbLRT pvalues 
#' @export
#'
pbLRT <- function(happi_out, 
                  firth = TRUE, 
                  spline_df = 3, 
                  max_iterations = 1000, 
                  change_threshold = 0.01, 
                  epsilon = 0, 
                  nstarts = 1, 
                  h0_param = 2, 
                  B = 1000){
 LL_alt <-  tail(happi_out$loglik$loglik[!is.na(happi_out$loglik$loglik)],1)
 LL_null <- tail(happi_out$loglik$loglik_null[!is.na(happi_out$loglik$loglik_null)],1)
 LL_alt_nopenalty <- tail(happi_out$loglik$loglik_nopenalty[!is.na(happi_out$loglik$loglik_nopenalty)],1)
 LL_null_nopenalty <- tail(happi_out$loglik$loglik_null_nopenalty[!is.na(happi_out$loglik$loglik_null_nopenalty)],1)
 
 beta_null <- as.vector(tail(happi_out$beta_null[!is.na(happi_out$beta_null[,1]),],1)[1])
 f_null <- as.vector(tail(happi_out$f_null[!is.na(happi_out$f_null[,1]),],1))
 quality_var <- happi_out$quality_var
 cbind(f_null, quality_var)
 covariate <- happi_out$covariate
 covariate_null <- covariate[, -h0_param]
 if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=length(covariate_null))
 
 t.observed <- 2 * (LL_alt - LL_null)
 t.observed_nopenalty <- 2 * (LL_alt_nopenalty-LL_null_nopenalty)    
 
   # lambda_props <- expit(covariate_null %*% beta_null)
   # colnames(lambda_props) <- c("lambda_props")
   # my_table <- cbind(f_null,quality_var,lambda_props) %>% as.data.frame()
   # # Draw from a Bernoulli distribution with p = lambda for both groups 
   # my_table$lambda <- rbinom(nrow(my_table), 1, my_table$lambda_props)
   # 
   # boot_gamma.table$Y_prop <- ifelse(boot_gamma.table$lambda == 1, 1-exp(-null_gamma1*boot_gamma.table$coverage), 0)
   # boot_gamma.table$Y <- rbinom(nrow(draw), 1, boot_gamma.table$Y_prop)
   # my_bootdata <- boot_gamma.table 
   # 
   # Y0 <- my_bootdata %>% dplyr::filter(Y == 0)
   # Y1 <- my_bootdata %>% dplyr::filter(Y == 1)
   # Mi0 = as.matrix(Y0['coverage'])
   # Mi0_mean = mean(Mi0)
   # 
   # Mi1 = as.matrix(Y1['coverage'])
   # Mi1_mean = mean(Mi1) 
   # 
   # gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
 
}
 
 


