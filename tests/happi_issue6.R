#### Look into whether activeSet double-dips or ignores weights 
# happi() has the weights = probs argument when updating f 
happi <- function(outcome,
                  covariate,
                  quality_var,
                  max_iterations = 50,
                  h0_param = 2,
                  nstarts = 1,
                  change_threshold = 0.05,
                  epsilon = 0,
                  method = "isotone",
                  spline_df = 4
) {
  
  # TODO(PT) take in formula
  
  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data
  
  nn <- length(outcome)
  
  stopifnot(nn == nrow(covariate) | nn == length(quality_var))
  
  pp <- ncol(covariate)
  
  if (ncol(covariate) > 2) warning("Amy hasn't properly checked that multiple covariates result in sensible output")
  
  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[my_order]
  outcome <- outcome[my_order]
  covariate <- covariate[my_order, ]
  
  covariate_null <- covariate[, -h0_param]
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=nn)
  
  calculate_p <- function(xbeta, ff) {
    pis <- ff * expit(xbeta) / (ff * expit(xbeta) + epsilon * (1 - expit(xbeta))) # for i: Y_i = 1
    pis_y0 <- (1 - ff) * expit(xbeta) / ((1 - ff) * expit(xbeta) + (1 - epsilon) * (1 - expit(xbeta)))
    pis[outcome == 0] <- pis_y0[outcome == 0]
    pis
  }
  
  update_beta <- function(probs, covariate) {
    logistf(probs ~ covariate - 1)$coef
  }
  
  
  update_f <- function(probs, tuning_param = 50,
                       method = "isotone",
                       spline_df = 4) {
    
    if(method == "isotone"){
      loss_fn <- function(x) -1 * sum(probs * (outcome * x - log(1 + exp(x)))) + sum(cosh((x / tuning_param)^2))
      loss_gradient <-  function(x) -1 * (probs * (outcome - exp(x) / (1 + exp(x)))) + (2 * x / tuning_param) * sinh((x / tuning_param)^2)
      
      ff_estimate <- activeSet(isomat = cbind(1:(nn-1), 2:nn), # define monotonicity
                               mySolver = fSolver,
                               fobj = loss_fn,
                               gobj = loss_gradient,
                               y = outcome,
                               weights = probs) ## this doesn't double dip on weights? (since they are specified in the fn and grad?)
      
      return(ff_estimate$x)
    } else if(method == "spline"){
      spline_basis <- cbind(1,iSpline(quality_var, df= spline_df, degree = 2, intercept = TRUE))
      b_start <- numeric(ncol(spline_basis))
      spline_criterion <- function(b){
        logit_means <- rowSums(do.call(cbind,lapply(1:length(b),
                                                    function(k) b[k]*spline_basis[,k,drop = FALSE])))
        return(-1*sum(probs*(outcome*logit_means - log(1 + exp(logit_means)))))
      }
      spline_fit <- optim(b_start,spline_criterion,method = "L-BFGS-B",
                          lower = c(-Inf,rep(0,length(b_start) - 1)),
                          upper = rep(Inf, length(b_start)))
      
      best_b <- spline_fit$par
      fitted_f_tilde <-
        rowSums(do.call(cbind,lapply(1:length(best_b),
                                     function(k)
                                       best_b[k]*
                                       spline_basis[,k,drop = FALSE])))
      return(fitted_f_tilde)
    }
  }
  
  
  incomplete_loglik <- function(xbeta, ff) {
    
    prob_lambda <- expit(xbeta)
    
    sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
               (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
  }
  
  ## no multiple starts for now
  my_estimates <- tibble("iteration" = 0:max_iterations,
                         "epsilon" = epsilon,
                         "loglik" = NA,
                         "loglik_null" = NA)
  
  my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
  my_estimated_beta_null <- matrix(NA, nrow = max_iterations + 1, ncol = pp - 1)
  my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_fitted_xbeta_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  
  my_estimated_beta[1, ] <- rep(0, pp)
  my_estimated_beta_null[1, ] <- rep(0, pp - 1)
  
  my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
  my_fitted_xbeta_null[1, ] <- c(covariate_null %*% my_estimated_beta_null[1, ])
  
  my_estimated_f[1, ] <- rep(mean(outcome), nn)
  my_estimated_f_null[1, ] <- rep(mean(outcome), nn)
  # f-tilde = logit(f)
  my_estimated_ftilde[1, ] <- logit(my_estimated_f[1, ])
  my_estimated_ftilde_null[1, ] <- logit(my_estimated_f_null[1, ])
  
  my_estimated_p[1, ] <- calculate_p(xbeta = my_fitted_xbeta[1, ],
                                     ff = my_estimated_f[1, ])
  my_estimated_p_null[1, ] <- calculate_p(xbeta = my_fitted_xbeta_null[1, ],
                                          ff = my_estimated_f_null[1, ])
  
  my_estimates[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                 ff = my_estimated_f[1, ])
  my_estimates[1, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                      ff = my_estimated_f_null[1, ])
  
  tt <- 1
  keep_going <- TRUE
  while (tt <= max_iterations & keep_going) {
    
    tt <- tt + 1
    
    ### alternative
    my_estimated_beta[tt, ] <- update_beta(probs=my_estimated_p[tt - 1, ], covariate=covariate)
    my_fitted_xbeta[tt, ] <- c(covariate %*% my_estimated_beta[tt, ])
    
    my_estimated_ftilde[tt, ] <- update_f(probs=my_estimated_p[tt - 1, ],
                                          method = method,
                                          spline_df = spline_df)
    my_estimated_f[tt, ] <- expit(my_estimated_ftilde[tt, ])
    
    my_estimated_p[tt, ] <- calculate_p(xbeta=my_fitted_xbeta[tt, ], ff=my_estimated_f[tt, ])
    
    my_estimates[tt, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[tt, ],
                                                    ff = my_estimated_f[tt, ])
    
    ### null
    my_estimated_beta_null[tt, ] <- update_beta(probs=my_estimated_p_null[tt - 1, ], covariate=covariate_null)
    my_fitted_xbeta_null[tt, ] <- c(covariate_null %*% my_estimated_beta_null[tt, ])
    
    my_estimated_ftilde_null[tt, ] <- update_f(probs=my_estimated_p_null[tt - 1, ],
                                               method = method,
                                               spline_df = spline_df)
    my_estimated_f_null[tt, ] <- expit(my_estimated_ftilde_null[tt, ])
    
    my_estimated_p_null[tt, ] <- calculate_p(xbeta=my_fitted_xbeta_null[tt, ], ff=my_estimated_f_null[tt, ])
    
    my_estimates[tt, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[tt, ],
                                                         ff = my_estimated_f_null[tt, ])
    
    ## maybe just log-likelihood changing?
    if (tt > 15) {
      # change_beta <- max(abs((my_estimated_beta[tt, ] - my_estimated_beta[tt - 1, ])/my_estimated_beta[tt - 1, ]))
      # change_beta_null <- max(abs((my_estimated_beta_null[tt, ] - my_estimated_beta_null[tt - 1, ])/my_estimated_beta_null[tt - 1, ]))
      # change_f <- max(abs((my_estimated_f[tt, ] - my_estimated_f[tt - 1, ])/my_estimated_f[tt - 1, ]))
      # change_f_null <- max(abs((my_estimated_f_null[tt, ] - my_estimated_f_null[tt - 1, ])/my_estimated_f_null[tt - 1, ]))
      
      pct_change_llks <- 100*max(abs((my_estimates[(tt - 4):tt, "loglik_null"] - my_estimates[(tt - 5):(tt - 1), "loglik_null"])/my_estimates[(tt - 5):(tt - 1), "loglik_null"]),
                                 abs((my_estimates[(tt - 4):tt, "loglik"] - my_estimates[(tt - 5):(tt - 1), "loglik"])/my_estimates[(tt - 5):(tt - 1), "loglik"]))
      keep_going <- pct_change_llks > change_threshold
      
      if(!keep_going) message(paste("Converged after", tt, "iterations; LL % change:", round(pct_change_llks, 3)))
    }
    
  }
  if (tt == max_iterations) {
    message(paste("Had not converged after", tt, "iterations; LL % change:", round(pct_change_llks, 3)))
  }
  
  my_estimates$LRT <- 2*(my_estimates$loglik - my_estimates$loglik_null)
  my_estimates$pvalue <- 1 - pchisq(my_estimates$LRT, df=1)
  
  return(list("loglik" = my_estimates,
              "beta" = my_estimated_beta,
              "beta_null" = my_estimated_beta_null,
              "f" = my_estimated_f,
              "f_null" = my_estimated_f_null,
              "p" = my_estimated_p,
              "p_null" = my_estimated_p_null,
              "quality_var" = quality_var,
              "outcome" = outcome,
              "covariate" = covariate))
  
  
}
## happi_mod() takes out the weights = probs argument when updating f 
happi_mod <- function(outcome,
                  covariate,
                  quality_var,
                  max_iterations = 50,
                  h0_param = 2,
                  nstarts = 1,
                  change_threshold = 0.05,
                  epsilon = 0,
                  method = "isotone",
                  spline_df = 4
) {
  
  # TODO(PT) take in formula
  
  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data
  
  nn <- length(outcome)
  
  stopifnot(nn == nrow(covariate) | nn == length(quality_var))
  
  pp <- ncol(covariate)
  
  if (ncol(covariate) > 2) warning("Amy hasn't properly checked that multiple covariates result in sensible output")
  
  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[my_order]
  outcome <- outcome[my_order]
  covariate <- covariate[my_order, ]
  
  covariate_null <- covariate[, -h0_param]
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=nn)
  
  calculate_p <- function(xbeta, ff) {
    pis <- ff * expit(xbeta) / (ff * expit(xbeta) + epsilon * (1 - expit(xbeta))) # for i: Y_i = 1
    pis_y0 <- (1 - ff) * expit(xbeta) / ((1 - ff) * expit(xbeta) + (1 - epsilon) * (1 - expit(xbeta)))
    pis[outcome == 0] <- pis_y0[outcome == 0]
    pis
  }
  
  update_beta <- function(probs, covariate) {
    logistf(probs ~ covariate - 1)$coef
  }
  
  
  update_f <- function(probs, tuning_param = 50,
                       method = "isotone",
                       spline_df = 4) {
    
    if(method == "isotone"){
      loss_fn <- function(x) -1 * sum(probs * (outcome * x - log(1 + exp(x)))) + sum(cosh((x / tuning_param)^2))
      loss_gradient <-  function(x) -1 * (probs * (outcome - exp(x) / (1 + exp(x)))) + (2 * x / tuning_param) * sinh((x / tuning_param)^2)
      
      ff_estimate <- activeSet(isomat = cbind(1:(nn-1), 2:nn), # define monotonicity
                               mySolver = fSolver,
                               fobj = loss_fn,
                               gobj = loss_gradient,
                               y = outcome)
      
      return(ff_estimate$x)
    } else if(method == "spline"){
      spline_basis <- cbind(1,iSpline(quality_var, df= spline_df, degree = 2, intercept = TRUE))
      b_start <- numeric(ncol(spline_basis))
      spline_criterion <- function(b){
        logit_means <- rowSums(do.call(cbind,lapply(1:length(b),
                                                    function(k) b[k]*spline_basis[,k,drop = FALSE])))
        return(-1*sum(probs*(outcome*logit_means - log(1 + exp(logit_means)))))
      }
      spline_fit <- optim(b_start,spline_criterion,method = "L-BFGS-B",
                          lower = c(-Inf,rep(0,length(b_start) - 1)),
                          upper = rep(Inf, length(b_start)))
      
      best_b <- spline_fit$par
      fitted_f_tilde <-
        rowSums(do.call(cbind,lapply(1:length(best_b),
                                     function(k)
                                       best_b[k]*
                                       spline_basis[,k,drop = FALSE])))
      return(fitted_f_tilde)
    }
  }
  
  
  incomplete_loglik <- function(xbeta, ff) {
    
    prob_lambda <- expit(xbeta)
    
    sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
               (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
      sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                ff[outcome == 1] * prob_lambda[outcome == 1]))
  }
  
  ## no multiple starts for now
  my_estimates <- tibble("iteration" = 0:max_iterations,
                         "epsilon" = epsilon,
                         "loglik" = NA,
                         "loglik_null" = NA)
  
  my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
  my_estimated_beta_null <- matrix(NA, nrow = max_iterations + 1, ncol = pp - 1)
  my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_fitted_xbeta_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  
  my_estimated_beta[1, ] <- rep(0, pp)
  my_estimated_beta_null[1, ] <- rep(0, pp - 1)
  
  my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
  my_fitted_xbeta_null[1, ] <- c(covariate_null %*% my_estimated_beta_null[1, ])
  
  my_estimated_f[1, ] <- rep(mean(outcome), nn)
  my_estimated_f_null[1, ] <- rep(mean(outcome), nn)
  # f-tilde = logit(f)
  my_estimated_ftilde[1, ] <- logit(my_estimated_f[1, ])
  my_estimated_ftilde_null[1, ] <- logit(my_estimated_f_null[1, ])
  
  my_estimated_p[1, ] <- calculate_p(xbeta = my_fitted_xbeta[1, ],
                                     ff = my_estimated_f[1, ])
  my_estimated_p_null[1, ] <- calculate_p(xbeta = my_fitted_xbeta_null[1, ],
                                          ff = my_estimated_f_null[1, ])
  
  my_estimates[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                 ff = my_estimated_f[1, ])
  my_estimates[1, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                      ff = my_estimated_f_null[1, ])
  
  tt <- 1
  keep_going <- TRUE
  while (tt <= max_iterations & keep_going) {
    
    tt <- tt + 1
    
    ### alternative
    my_estimated_beta[tt, ] <- update_beta(probs=my_estimated_p[tt - 1, ], covariate=covariate)
    my_fitted_xbeta[tt, ] <- c(covariate %*% my_estimated_beta[tt, ])
    
    my_estimated_ftilde[tt, ] <- update_f(probs=my_estimated_p[tt - 1, ],
                                          method = method,
                                          spline_df = spline_df)
    my_estimated_f[tt, ] <- expit(my_estimated_ftilde[tt, ])
    
    my_estimated_p[tt, ] <- calculate_p(xbeta=my_fitted_xbeta[tt, ], ff=my_estimated_f[tt, ])
    
    my_estimates[tt, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[tt, ],
                                                    ff = my_estimated_f[tt, ])
    
    ### null
    my_estimated_beta_null[tt, ] <- update_beta(probs=my_estimated_p_null[tt - 1, ], covariate=covariate_null)
    my_fitted_xbeta_null[tt, ] <- c(covariate_null %*% my_estimated_beta_null[tt, ])
    
    my_estimated_ftilde_null[tt, ] <- update_f(probs=my_estimated_p_null[tt - 1, ],
                                               method = method,
                                               spline_df = spline_df)
    my_estimated_f_null[tt, ] <- expit(my_estimated_ftilde_null[tt, ])
    
    my_estimated_p_null[tt, ] <- calculate_p(xbeta=my_fitted_xbeta_null[tt, ], ff=my_estimated_f_null[tt, ])
    
    my_estimates[tt, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[tt, ],
                                                         ff = my_estimated_f_null[tt, ])
    
    ## maybe just log-likelihood changing?
    if (tt > 15) {
      # change_beta <- max(abs((my_estimated_beta[tt, ] - my_estimated_beta[tt - 1, ])/my_estimated_beta[tt - 1, ]))
      # change_beta_null <- max(abs((my_estimated_beta_null[tt, ] - my_estimated_beta_null[tt - 1, ])/my_estimated_beta_null[tt - 1, ]))
      # change_f <- max(abs((my_estimated_f[tt, ] - my_estimated_f[tt - 1, ])/my_estimated_f[tt - 1, ]))
      # change_f_null <- max(abs((my_estimated_f_null[tt, ] - my_estimated_f_null[tt - 1, ])/my_estimated_f_null[tt - 1, ]))
      
      pct_change_llks <- 100*max(abs((my_estimates[(tt - 4):tt, "loglik_null"] - my_estimates[(tt - 5):(tt - 1), "loglik_null"])/my_estimates[(tt - 5):(tt - 1), "loglik_null"]),
                                 abs((my_estimates[(tt - 4):tt, "loglik"] - my_estimates[(tt - 5):(tt - 1), "loglik"])/my_estimates[(tt - 5):(tt - 1), "loglik"]))
      keep_going <- pct_change_llks > change_threshold
      
      if(!keep_going) message(paste("Converged after", tt, "iterations; LL % change:", round(pct_change_llks, 3)))
    }
    
  }
  if (tt == max_iterations) {
    message(paste("Had not converged after", tt, "iterations; LL % change:", round(pct_change_llks, 3)))
  }
  
  my_estimates$LRT <- 2*(my_estimates$loglik - my_estimates$loglik_null)
  my_estimates$pvalue <- 1 - pchisq(my_estimates$LRT, df=1)
  
  return(list("loglik" = my_estimates,
              "beta" = my_estimated_beta,
              "beta_null" = my_estimated_beta_null,
              "f" = my_estimated_f,
              "f_null" = my_estimated_f_null,
              "p" = my_estimated_p,
              "p_null" = my_estimated_p_null,
              "quality_var" = quality_var,
              "outcome" = outcome,
              "covariate" = covariate))
  
  
}

expit <- function(input) {
  exp(input) / (1 + exp(input))
}

logit <- function(input) {
  log(input / (1 - input))
}


### example dataset 1
r27 <- structure(list(tongue = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                      mean_coverage = c(7.998170265,
                                        5.988207055, 14.94151039, 15.08557254, 4.050932633, 20.98467875,
                                        17.24238886, 1.422470366, 1.864070686, 6.578526287, 3.429256021,
                                        1.218447863, 20.17103072, 5.97223494, 9.673829539, 15.94390898,
                                        3.2883879, 3.13199016, 2.795140254, 2.29811946, 7.848127938,
                                        15.6438722, 4.197064174, 7.649882259, 3.153925241, 2.870245578,
                                        2.921719021, 4.982794781, 5.680954521, 1.492299764, 12.65457264,
                                        1.071999059, 2.557173838, 1.847580956, 5.456339015, 12.03145172,
                                        18.40454984, 2.829514291, 26.3464049, 7.841277141, 2.448978134,
                                        14.69551371, 13.09895311),
                      `Ribosomal protein L27` = c(1, 1,
                                                  1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                  1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1)),
                 class = "data.frame", row.names = c(NA,
                                                     -43L))


## Let's do a first pass comparison when we change the activeSet function to see if it double dips/produces different results 
happi_r27 <- happi(outcome = r27$`Ribosomal protein L27`,
                   covariate = cbind(1, r27$tongue),
                   quality_var = r27$mean_coverage,
                   max_iterations = 500,
                   nstarts = 1,
                   epsilon = 0) # this is when weights = probs
happi_r27_mod <- happi_mod(outcome = r27$`Ribosomal protein L27`,
                           covariate = cbind(1, r27$tongue),
                           quality_var = r27$mean_coverage,
                           max_iterations = 500,
                           nstarts = 1,
                           epsilon = 0) # this is when we take out weights = probs
# okay let's compare results with both approaches to updating f 
happi_r27$beta %>% tail
happi_r27_mod$beta %>% tail # looks the same

happi_r27$loglik %>% tail(1) 
happi_r27_mod$loglik %>% tail(1) # identical results 

happi_r27$f %>% tail(2) %>% logit 
happi_r27_mod$f %>% tail(2) %>% logit # identical results
happi_r27$f %>% tail(1) %>% logit %>% rbind(happi_r27_mod$f %>% tail(1) %>% logit) # matched

happi_r27$f %>% head(2) %>% logit
happi_r27_mod$f %>% head(2) %>% logit # identical results 

happi_r27$p %>% tail()
happi_r27_mod$p %>% tail() # identical results

happi_r27$loglik$pvalue %>% tail()
happi_r27_mod$loglik$pvalue %>% tail() # identical results


### example dataset 2 
mem <- structure(list(tongue = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                      mean_coverage = c(7.998170265,
                                        5.988207055, 14.94151039, 15.08557254, 4.050932633, 20.98467875,
                                        17.24238886, 1.422470366, 1.864070686, 6.578526287, 3.429256021,
                                        1.218447863, 20.17103072, 5.97223494, 9.673829539, 15.94390898,
                                        3.2883879, 3.13199016, 2.795140254, 2.29811946, 7.848127938,
                                        15.6438722, 4.197064174, 7.649882259, 3.153925241, 2.870245578,
                                        2.921719021, 4.982794781, 5.680954521, 1.492299764, 12.65457264,
                                        1.071999059, 2.557173838, 1.847580956, 5.456339015, 12.03145172,
                                        18.40454984, 2.829514291, 26.3464049, 7.841277141, 2.448978134,
                                        14.69551371, 13.09895311),
                      `Membrane protein insertase Oxa1/YidC/SpoIIIJ, required for the localization of integral membrane proteins` =
                        c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1)),
                 class = "data.frame", row.names = c(NA, -43L))


names(mem)[3] <- "menbraneprotein"
happi_mem <- happi(outcome = mem$menbraneprotein,
                   covariate = cbind(1, mem$tongue),
                   quality_var = mem$mean_coverage,
                   max_iterations = 100,
                   nstarts = 1,
                   epsilon = 0)
happi_mem_mod <- happi(outcome = mem$menbraneprotein,
                   covariate = cbind(1, mem$tongue),
                   quality_var = mem$mean_coverage,
                   max_iterations = 100,
                   nstarts = 1,
                   epsilon = 0)

happi_mem$beta %>% tail 
happi_mem_mod$beta %>% tail # match

happi_mem$f %>% tail(1) %>% logit %>% rbind(happi_mem_mod$f %>% tail(1) %>% logit) # match

happi_mem$loglik %>% tail 
happi_mem_mod$loglik %>% tail # match

happi_mem$loglik$loglik %>% tail(20) 
happi_mem_mod$loglik$loglik %>% tail(20)

happi_mem$loglik$loglik_null %>% tail(20) 
happi_mem_mod$loglik$loglik_null %>% tail(20) 

happi_mem$loglik$pvalue %>% tail(20)
happi_mem_mod$loglik$pvalue %>% tail(20)

## Based on these two it looks like whether we specify the weights or not doesn't produce different results of f and that there isn't any impact on the other 
## estimated parameters/results 