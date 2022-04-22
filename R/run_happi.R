#' Main function for p=q=1
#'
#' @param outcome length-n vector
#' @param covariate n x p matrix
#' @param quality_var length-n vector  TODO(turn into n x q matrix)
#' @param max_iterations the maximum number of EM steps that the algorithm will run for
#' @param min_iterations the minimum number of EM steps that the algorithm will run for
#' @param h0_param the column index in covariate that has beta=zero under the null
#' @param nstarts number of starts TODO
#' @param change_threshold algorithm will terminate early if the likelihood changes by this percentage or less for 5 iterations in a row for both the alternative and the null
#' @param epsilon prob of ...
#' @param method method for estimating f. Defaults to "isotone" for isotonic
#' regression fit; "spline" fits a monotone spline with df determined by
#' argument spline_df
#' @param random_starts whether to pick the starting values of beta's randomly. Defaults to FALSE.
#' @param firth use firth penalty? Default is TRUE.
#' @param spline_df degrees of freedom (in addition to intercept) to use in
#' monotone spline fit
#'
#'
#' @import tibble
#' @import stats
#' @import splines2
#' @importFrom isotone activeSet fSolver
#' @importFrom logistf logistf
#'
#' @export
happi <- function(outcome,
                  covariate,
                  quality_var,
                  max_iterations = 50,
                  min_iterations = 15,
                  h0_param = 2,
                  nstarts = 1,
                  change_threshold = 0.05,
                  epsilon = 0,
                  method = "isotone",
                  random_starts = FALSE,
                  firth = TRUE,
                  spline_df = 4
) {


  # TODO(PT) take in formula

  stopifnot(all(!is.na(c(outcome, covariate, quality_var)))) # some missing data

  nn <- length(outcome)

  stopifnot(nn == nrow(covariate) | nn == length(quality_var))

  pp <- ncol(covariate)

  if (ncol(covariate) > 2) warning("Amy hasn't properly checked that multiple covariates result in sensible output")
  if(h0_param != 2) warning("Amy hasn't properly checked that testing a different parameter results in sensible output")

  ## reorder all elements of all data by ordering in quality_var
  ## TODO(change back at end)
  my_order <- order(quality_var)
  quality_var <- quality_var[my_order]
  outcome <- outcome[my_order]
  covariate <- covariate[my_order, ]

  if (pp == 1) {
    covariate <- matrix(covariate, ncol = 1, nrow = nn)
    covariate_null <- matrix(1, ncol = 1, nrow = nn)
  } else {
    covariate_null <- covariate[, -h0_param]
  }
  if (!is.matrix(covariate_null)) covariate_null <- matrix(covariate_null, nrow=nn)

  calculate_p <- function(xbeta, ff) {
    pis <- ff * expit(xbeta) / (ff * expit(xbeta) + epsilon * (1 - expit(xbeta))) # for i: Y_i = 1
    pis_y0 <- (1 - ff) * expit(xbeta) / ((1 - ff) * expit(xbeta) + (1 - epsilon) * (1 - expit(xbeta)))
    pis[outcome == 0] <- pis_y0[outcome == 0]
    pis
  }

  update_beta <- function(probs, covariate) {

    if (!firth) {
      # glm(probs ~ covariate - 1, family=binomial)$coef

      ## prevents warnings about `non-integer #successes in a binomial glm!`
      ## doesn't alter coefficient estimates compared to "binomial", only std errors, which we don't use
      coefs <- glm(probs ~ covariate - 1, family= quasibinomial)$coef
    } else {
      coefs <- logistf(probs ~ covariate - 1)$coef
    }
    coefs
  }


  update_f <- function(probs,
                       tuning_param = 50,
                       method = "isotone",
                       spline_df = 4) {

    if(method == "isotone") {
      loss_fn <- function(x) -1 * sum(probs * (outcome * x - log(1 + exp(x)))) + sum(cosh((x / tuning_param)^2))
      loss_gradient <-  function(x) -1 * (probs * (outcome - exp(x) / (1 + exp(x)))) + (2 * x / tuning_param) * sinh((x / tuning_param)^2)

      ff_estimate <- activeSet(isomat = cbind(1:(nn-1), 2:nn), # define monotonicity
                               mySolver = fSolver,
                               fobj = loss_fn,
                               gobj = loss_gradient)
      return(ff_estimate$x)
    } else if (method == "spline") {

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
      fitted_f_tilde <- rowSums(do.call(cbind,lapply(1:length(best_b),
                                                     function(k)
                                                       best_b[k]*spline_basis[,k,drop = FALSE])))
      return(fitted_f_tilde)
    } else {
      stop("Invalid input to `method`. Choose 'isotone' or 'spline'.")
    }
  }


  incomplete_loglik <- function(xbeta, ff, firth = TRUE) {

    prob_lambda <- expit(xbeta)

    if (!firth) {
      sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
                 (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
        sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                  ff[outcome == 1] * prob_lambda[outcome == 1]))
    } else {
      ll <-  sum(log( (1 - epsilon) * (1 - prob_lambda[outcome == 0]) +
                        (1 - ff[outcome == 0]) * prob_lambda[outcome == 0])) +
        sum(log(epsilon * (1 - prob_lambda[outcome == 1]) +
                  ff[outcome == 1] * prob_lambda[outcome == 1]))
      penalty <- 0.5*det(t(covariate)%*%diag(as.numeric(prob_lambda)*(
        1 - as.numeric(prob_lambda)))%*%covariate)
      # print(c(ll, penalty))
      return(ll - penalty)
    }
  }

  ## no multiple starts for now
  my_estimates <- tibble("iteration" = 0:max_iterations,
                         "epsilon" = epsilon,
                         "loglik" = NA,
                         "loglik_null" = NA)

  my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
  my_estimated_beta_null <- matrix(NA, nrow = max_iterations + 1, ncol = max(1, pp - 1))
  my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_fitted_xbeta_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_f_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_ftilde_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
  my_estimated_p_null <- matrix(NA, nrow = max_iterations + 1, ncol = nn)

  if (random_starts) {
    my_estimated_beta[1, ] <- rnorm(pp)
    my_estimated_beta_null[1, ] <- rnorm(pp - 1)
  } else {
    my_estimated_beta[1, ] <- rep(0, pp)
    my_estimated_beta_null[1, ] <- rep(0, max(1, pp - 1))
  }

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
                                                 ff = my_estimated_f[1, ],
                                                 firth = firth)
  my_estimates[1, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[1, ],
                                                      ff = my_estimated_f_null[1, ],
                                                      firth = firth)

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
                                                    ff = my_estimated_f[tt, ],
                                                    firth = firth)

    ### null
    my_estimated_beta_null[tt, ] <- update_beta(probs=my_estimated_p_null[tt - 1, ], covariate=covariate_null)
    my_fitted_xbeta_null[tt, ] <- c(covariate_null %*% my_estimated_beta_null[tt, ])

    my_estimated_ftilde_null[tt, ] <- update_f(probs=my_estimated_p_null[tt - 1, ],
                                               method = method,
                                               spline_df = spline_df)
    my_estimated_f_null[tt, ] <- expit(my_estimated_ftilde_null[tt, ])

    my_estimated_p_null[tt, ] <- calculate_p(xbeta=my_fitted_xbeta_null[tt, ], ff=my_estimated_f_null[tt, ])

    my_estimates[tt, "loglik_null"] <- incomplete_loglik(xbeta = my_fitted_xbeta_null[tt, ],
                                                         ff = my_estimated_f_null[tt, ],
                                                         firth = firth)

    ## maybe just log-likelihood changing?
    if ((tt > min_iterations) & (my_estimates[tt, "loglik"] > my_estimates[tt, "loglik_null"])) {

      pct_change_llks <- 100*max(abs((my_estimates[(tt - 4):tt, "loglik_null"] - my_estimates[(tt - 5):(tt - 1), "loglik_null"])/my_estimates[(tt - 5):(tt - 1), "loglik_null"]),
                                 abs((my_estimates[(tt - 4):tt, "loglik"] - my_estimates[(tt - 5):(tt - 1), "loglik"])/my_estimates[(tt - 5):(tt - 1), "loglik"]))
      keep_going <- pct_change_llks > change_threshold

      if(!keep_going) message(paste("Converged after", tt, "iterations; LL % change:", round(pct_change_llks, 3)))
    }

  }

  tt_restart <- 1
  if (my_estimates[tt, "loglik"] < my_estimates[tt, "loglik_null"]) {

    message("Likelihood greater under null; restarting...")
    my_estimated_beta <- matrix(NA, nrow = max_iterations + 1, ncol = pp)
    my_fitted_xbeta <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_f <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_ftilde <- matrix(NA, nrow = max_iterations + 1, ncol = nn)
    my_estimated_p <- matrix(NA, nrow = max_iterations + 1, ncol = nn)

    if (pp > 1) {
      my_estimated_beta[1, 1] <- my_estimated_beta_null[tt] ## start at converged null
      my_estimated_beta[1, h0_param] <- 0
    } else {
      my_estimated_beta[1, 1] <- 0
    }
    stopifnot(h0_param == 2)
    my_fitted_xbeta[1, ] <- c(covariate %*% my_estimated_beta[1, ])
    my_estimated_f[1, ] <- my_estimated_f_null[tt, ]
    my_estimated_ftilde[1, ] <- logit(my_estimated_f[1, ])
    my_estimated_p[1, ] <- calculate_p(xbeta = my_fitted_xbeta[1, ],
                                       ff = my_estimated_f[1, ])
    my_estimates[1, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[1, ],
                                                   ff = my_estimated_f[1, ])


    ## restart at null model if likelihood greater under null than alternative

    keep_going <- TRUE
    while (tt_restart <= max_iterations & keep_going) {

      tt_restart <- tt_restart + 1

      ### alternative
      my_estimated_beta[tt_restart, ] <- update_beta(probs=my_estimated_p[tt_restart - 1, ], covariate=covariate)
      my_fitted_xbeta[tt_restart, ] <- c(covariate %*% my_estimated_beta[tt_restart, ])

      my_estimated_ftilde[tt_restart, ] <- update_f(probs=my_estimated_p[tt_restart - 1, ],
                                                    method = method,
                                                    spline_df = spline_df)
      my_estimated_f[tt_restart, ] <- expit(my_estimated_ftilde[tt_restart, ])

      my_estimated_p[tt_restart, ] <- calculate_p(xbeta=my_fitted_xbeta[tt_restart, ], ff=my_estimated_f[tt_restart, ])

      my_estimates[tt_restart, "loglik"] <- incomplete_loglik(xbeta = my_fitted_xbeta[tt_restart, ],
                                                              ff = my_estimated_f[tt_restart, ])

      if ((tt_restart > min_iterations) & (my_estimates[tt_restart, "loglik"] > my_estimates[tt, "loglik_null"])) {

        pct_change_llks <- 100*abs((my_estimates[(tt_restart - 4):tt_restart, "loglik"] - my_estimates[(tt_restart - 5):(tt_restart - 1), "loglik"])/my_estimates[(tt_restart - 5):(tt_restart - 1), "loglik"])
        keep_going <- pct_change_llks > change_threshold

        if(!keep_going) message(paste("Converged after", tt_restart, "iterations; LL % change:", round(pct_change_llks, 3)))
      }

    }

  }

  if (tt_restart == max_iterations + 1) {
    message("Restarting to estimate beta_alt didn't work. Not sure what's happening...")
    # message(paste("Had not converged after", tt_restart - 1, "iterations; LL % change:", round(pct_change_llks, 3)))
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

