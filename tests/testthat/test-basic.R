context("basic checks that the package works")
library(happi)

test_that("happi runs", {

  set.seed(2)
  nn <- 20
  mm <- seq(30, 100, length.out = nn)
  xx <- cbind(1, rep(c(0, 1), nn/2))
  beta <- c(0, 0)
  epsilon <- 0

  true_lambda_probs <- expit(c(xx %*% beta))
  lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)

  ## fake f: expit(gamma0 + gamma1 * m)
  ## goal: Pr(Y = 1 | lambda = 1, m = 30) = 0.2 = expit(gamma0 + gamma1 * 30)
  ## goal: Pr(Y = 1 | lambda = 1, m = 100) = 0.8 = expit(gamma0 + gamma1 * 100)
  logit(0.2) # -1.38
  logit(0.8) # 1.38
  gamma1 <- (logit(0.8) - logit(0.2)) / (100 - 30)
  gamma0 <- logit(0.8) - gamma1 * 100

  true_f <- expit(gamma0 + gamma1 * mm)

  true_prob_y_equal_1 <- true_f # correct for lambdas == 1
  true_prob_y_equal_1[lambdas == 0] <- epsilon

  ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

  devtools::load_all()
  hh0 <- happi(outcome = ys,
              covariate = xx,
              quality_var = mm,
              max_iterations = 30,
              nstarts = 1,
              epsilon = 0)
  hh0$loglik$loglik
  hh0$loglik$loglik_null
  # 1e5 * (hh0$loglik$loglik %>% exp)
  # 1e5 * (hh5$loglik$loglik %>% exp)
  # hh5$loglik$loglik %>% plot

  hh5 <- happi(outcome = ys,
          covariate = xx,
          quality_var = mm,
          max_iterations = 10,
          nstarts = 1,
          epsilon = 0.05)



})
