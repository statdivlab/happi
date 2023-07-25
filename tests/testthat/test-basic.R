
test_that("happi runs", {

  set.seed(21)
  nn <- 100
  mm <- seq(30, 100, length.out = nn)
  xx <- cbind(1, rep(c(0, 1), nn/2))
  beta <- c(0, 0)
  epsilon <- 0

  true_lambda_probs <- happi::expit(c(xx %*% beta))
  lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)

  ## fake f: expit(gamma0 + gamma1 * m)
  ## goal: Pr(Y = 1 | lambda = 1, m = 30) = 0.2 = expit(gamma0 + gamma1 * 30)
  ## goal: Pr(Y = 1 | lambda = 1, m = 100) = 0.8 = expit(gamma0 + gamma1 * 100)
  logit(0.2) # -1.38
  logit(0.8) # 1.38
  gamma1 <- (logit(0.8) - logit(0.2)) / (100 - 30)
  gamma0 <- logit(0.8) - gamma1 * 100

  true_f <- happi::expit(gamma0 + gamma1 * mm)

  true_prob_y_equal_1 <- true_f # correct for lambdas == 1
  true_prob_y_equal_1[lambdas == 0] <- epsilon

  ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

  hh0_isotone <- happi::happi(outcome = ys,
                       covariate = xx,
                       quality_var = mm,
                       max_iterations = 1000,
                       method = "isotone",
                       nstarts = 1,
                       epsilon = 0, 
                       change_threshold=0.1)

  hh0_spline <- happi::happi(outcome = ys,
                       covariate = xx,
                       quality_var = mm,
                       max_iterations = 1000,
                       method = "splines",
                       nstarts = 1,
                       epsilon = 0)
  
  # check that the epsilon argument will take in either a single value or a vector of length n 
  hh0_spline_mult_eps <- happi::happi(outcome = ys,
                                      covariate = xx,
                                      quality_var = mm,
                                      max_iterations = 1000,
                                      method = "splines",
                                      nstarts = 1,
                                      epsilon = rep(0, length(ys)))

  expect_type(hh0_isotone, "list")
  expect_type(hh0_spline, "list")
  expect_equal(hh0_spline$loglik$pvalue, hh0_spline_mult_eps$loglik$pvalue)

})


test_that("happi works with formulas", {
  
  set.seed(21)
  nn <- 100
  mm <- seq(30, 100, length.out = nn)
  xx <- cbind(1, rep(c(0, 1), nn/2))
  beta <- c(0, 0)
  epsilon <- 0
  
  true_lambda_probs <- happi::expit(c(xx %*% beta))
  lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
  
  ## fake f: expit(gamma0 + gamma1 * m)
  ## goal: Pr(Y = 1 | lambda = 1, m = 30) = 0.2 = expit(gamma0 + gamma1 * 30)
  ## goal: Pr(Y = 1 | lambda = 1, m = 100) = 0.8 = expit(gamma0 + gamma1 * 100)
  logit(0.2) # -1.38
  logit(0.8) # 1.38
  gamma1 <- (logit(0.8) - logit(0.2)) / (100 - 30)
  gamma0 <- logit(0.8) - gamma1 * 100
  
  true_f <- happi::expit(gamma0 + gamma1 * mm)
  
  true_prob_y_equal_1 <- true_f # correct for lambdas == 1
  true_prob_y_equal_1[lambdas == 0] <- epsilon
  
  ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)
  
  res <- happi::happi(outcome = ys,
                      covariate = xx,
                      quality_var = mm,
                      max_iterations = 1000,
                      method = "splines",
                      nstarts = 1,
                      epsilon = 0)
  
  dat <- data.frame(cov = xx[, 2],
                    qual_var = mm)
  res_form <- happi::happi(outcome = ys,
                           covariate_formula = ~ cov,
                           covariate_formula_h0 = ~ 1,
                           quality_var_formula = ~ qual_var,
                           data = dat,
                           max_iterations = 1000,
                           method = "splines",
                           nstarts = 1,
                           epsilon = 0)
  
  # check that both forms of data give same results
  expect_true(all.equal(res$loglik, res_form$loglik))
})

# test_that("happi works with multiple covariates", {
#   
#   set.seed(21)
#   nn <- 100
#   mm <- seq(30, 100, length.out = nn)
#   xx <- cbind(1, rep(c(0, 1), nn/2))
#   beta <- c(0, 0)
#   epsilon <- 0
#   
#   true_lambda_probs <- happi::expit(c(xx %*% beta))
#   lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
#   
#   ## fake f: expit(gamma0 + gamma1 * m)
#   ## goal: Pr(Y = 1 | lambda = 1, m = 30) = 0.2 = expit(gamma0 + gamma1 * 30)
#   ## goal: Pr(Y = 1 | lambda = 1, m = 100) = 0.8 = expit(gamma0 + gamma1 * 100)
#   logit(0.2) # -1.38
#   logit(0.8) # 1.38
#   gamma1 <- (logit(0.8) - logit(0.2)) / (100 - 30)
#   gamma0 <- logit(0.8) - gamma1 * 100
#   
#   true_f <- happi::expit(gamma0 + gamma1 * mm)
#   
#   true_prob_y_equal_1 <- true_f # correct for lambdas == 1
#   true_prob_y_equal_1[lambdas == 0] <- epsilon
#   
#   ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)
#   
#   dat <- data.frame(cov = xx[, 2],
#                     qual_var = mm,
#                     rand1 = rnorm(nn),
#                     rand2 = rnorm(nn))
#   
#   res <- happi::happi_multi_cov(outcome = ys,
#                            covariate_formula = ~ cov + rand1,
#                            covariate_formula_h0 = ~ rand1,
#                            quality_var_formula = ~ qual_var,
#                            data = dat,
#                            max_iterations = 1000,
#                            method = "splines",
#                            nstarts = 1,
#                            epsilon = 0)
#   
#   # check that both forms of data give same results
#   expect_true(all.equal(res$loglik, res_form$loglik))
# })
