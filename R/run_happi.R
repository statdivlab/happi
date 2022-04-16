#' Main function for p=q=1
#'
#' @param outcome length-n vector
#' @param covariate n x p matrix
#' @param quality_var length-n vector  TODO(turn into n x q matrix)
#' @param max_iterations
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param nstarts number of starts TODO
#' @param epsilon prob of ...
#'
#' @import tibble
#' @importFrom isotone activeSet fSolver
#' @importFrom logistf logistf
#' @importFrom stats pchisq
#'
#' @export
happi <- function(outcome, covariate, quality_var,
                  max_iterations = 50,
                  h0_param = 2,
                  nstarts = 1,
                  epsilon = 0
) {

  # TODO(PT) take in formula

  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data

  nn <- length(outcome)

  stopifnot(nn == nrow(covariate) | nn == length(quality_var))

  pp <- ncol(covariate)

  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[order(quality_var)]
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


  update_f <- function(probs, tuning_param = 50) {

    loss_fn <- function(x) -1 * sum(probs * (outcome * x - log(1 + exp(x)))) + sum(cosh((x / tuning_param)^2))
    loss_gradient <-  function(x) -1 * (probs * (outcome - exp(x) / (1 + exp(x)))) + (2 * x / tuning_param) * sinh((x / tuning_param)^2)

    ff_estimate <- activeSet(isomat = cbind(1:(nn-1), 2:nn), # define monotonicity
                             mySolver = fSolver,
                             fobj = loss_fn,
                             gobj = loss_gradient,
                             y = outcome,
                             weights = probs) ## this doesn't double dip on weights? (since they are specified in the fn and grad?)

    return(ff_estimate$x)
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

  my_estimated_beta[1, ] <- c(4, -2) # rep(0, pp)
  my_estimated_beta_null[1, ] <- c(4) #rep(0, pp - 1)

  my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
  my_fitted_xbeta_null[1, ] <- c(covariate_null %*% my_estimated_beta_null[1, ])

  my_estimated_f[1, ] <- rep(0.73, nn)
  my_estimated_f_null[1, ] <- rep(0.73, nn)
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
  other_convergence_conditions <- TRUE # TODO
  while (tt <= max_iterations & other_convergence_conditions) {

    tt <- tt + 1

    ### alternative
    my_estimated_beta[tt, ] <- update_beta(probs=my_estimated_p[tt - 1, ], covariate=covariate)
    my_fitted_xbeta[tt, ] <- c(covariate %*% my_estimated_beta[tt, ])

    my_estimated_ftilde[tt, ] <- update_f(probs=my_estimated_p[tt - 1, ])
    my_estimated_f[tt, ] <- expit(my_estimated_ftilde[tt, ])

    my_estimated_p[tt, ] <- calculate_p(xbeta=my_fitted_xbeta[tt, ], ff=my_estimated_f[tt, ])

    my_estimates[tt, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[tt, ],
                                                    ff = my_estimated_f[tt, ])

    ### null
    my_estimated_beta_null[tt, ] <- update_beta(probs=my_estimated_p_null[tt - 1, ], covariate=covariate_null)
    my_fitted_xbeta_null[tt, ] <- c(covariate_null %*% my_estimated_beta_null[tt, ])

    my_estimated_ftilde_null[tt, ] <- update_f(probs=my_estimated_p_null[tt - 1, ])
    my_estimated_f_null[tt, ] <- expit(my_estimated_ftilde_null[tt, ])

    my_estimated_p_null[tt, ] <- calculate_p(xbeta=my_fitted_xbeta_null[tt, ], ff=my_estimated_f_null[tt, ])

    my_estimates[tt, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[tt, ],
                                                         ff = my_estimated_f_null[tt, ])

  }

  my_estimates$LRT <- 2*(my_estimates$loglik - my_estimates$loglik_null)
  my_estimates$pvalue <- 1 - pchisq(my_estimates$LRT, df=1)

  return(list("loglik" = my_estimates,
              "beta" = my_estimated_beta,
              "beta_null" = my_estimated_beta_null,
              "f" = my_estimated_f,
              "f_null" = my_estimated_f_null,
              "p" = my_estimated_p,
              "p_null" = my_estimated_p_null))


}

