#' Function for checking update_f 
#'
#' @param probs outcomes p; probability between 0 and 1 
#' @param tuning_param tuning parameter. default set to  50 
#' @param method method for updating f. Can specify either  isotone or splines. default splines. 
#' @param spline_df degrees of freedom for splines. default 3. 
#' @param quality_var length-n vector; this is the quality variable vector, currently p = 1
#' @param nn length of outcome vector 
#' @param outcome length-n vector; this is the vector of a target gene's presence/absence; should be coded as 0 or 1 
#' @param tt E-M iteration value 
#'
#' @return warning for checking f update
#'
#' @examples
#' data("TM7_data")
#' # create design matrix
#' covariate <- model.matrix(~tongue, data = TM7_data) 
#' outcome <- TM7_data$`Ribosomal protein L27` 
#' quality_var <- TM7_data$mean_coverage
#' spline_df <- 3 
#' nn <- length(outcome)
#' probs <- rep(0.5,nn)
#' 
#' my_warningcheck_update_f <- warningcheck_update_f(probs = probs, 
#' tuning_param = 50, 
#' method = "splines", 
#' spline_df = spline_df, 
#' outcome = outcome, 
#' nn = nn, 
#' quality_var = quality_var)
#' @export
warningcheck_update_f <- function(probs,
                                  tuning_param = 50,
                                  method = "splines",
                                  spline_df = spline_df, 
                                  quality_var = quality_var, 
                                  nn = nn, 
                                  outcome = outcome, 
                                  tt = tt) {
  
  if(method == "isotone") {
    loss_fn <- function(x) -1 * sum(probs * (outcome * x - log(1 + exp(x)))) + sum(cosh((x / tuning_param)^2))
    loss_gradient <-  function(x) -1 * (probs * (outcome - exp(x) / (1 + exp(x)))) + (2 * x / tuning_param) * sinh((x / tuning_param)^2)
    
    ff_estimate <- activeSet(isomat = cbind(1:(nn-1), 2:nn), # define monotonicity
                             mySolver = fSolver,
                             fobj = loss_fn,
                             gobj = loss_gradient)
    return(list("fitted_f_tilde" = ff_estimate$x, "basis_weights" = NA))
  } else if (method %in% c("splines", "spline")) {
    
    spline_basis <- cbind(1, iSpline(quality_var, df= spline_df, degree = 2, intercept = TRUE))
    
    b_start <- rep(0, ncol(spline_basis))
    
    spline_criterion <- function(b) {
      logit_means <- rowSums(do.call(cbind, lapply(1:length(b),
                                                   function(k) b[k]*spline_basis[,k,drop = FALSE])))
      return(-1*sum(probs*(outcome*logit_means - log(1 + exp(logit_means)))))
    }
    
    spline_fit <- optim(b_start,
                        spline_criterion,
                        method = "L-BFGS-B",
                        # lower = c(-Inf, rep(0,length(b_start) - 1)),
                        # upper = rep(Inf, length(b_start))
                        lower = c(-1e7, rep(0,length(b_start) - 1)), ## to improve stability
                        upper = rep(100, length(b_start)) ## to improve stability
    )
    
    best_b <- spline_fit$par
    # see GitHub Issue #13 https://github.com/statdivlab/happi/issues/13
    if (any(abs(best_b - 100) < 0.01)) warning("spline basis weights close to boundaries; try reducing spline_df; warning at iteration",tt)
    if (any(abs(best_b > 1e4))) warning("spline basis weights close to boundaries; try reducing spline_df; warning at iteration",tt)
    fitted_f_tilde <- rowSums(do.call(cbind,lapply(1:length(best_b),
                                                   function(k)
                                                     best_b[k]*spline_basis[,k,drop = FALSE])))
    return(list("fitted_f_tilde" = fitted_f_tilde, "basis_weights" = best_b))
  } else {
    stop("Invalid input to `method`. Choose 'isotone' or 'spline'.")
  }
}