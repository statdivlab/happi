library(happi)
context("Test raotest")

test_that("happi fails to converge when expit(XB) is too large because optim needs finite values of 'fn'", {

  set.seed(21)
  nn <- 16
  mm <- seq(30, 100, length.out = nn)
  xx <- cbind(1, rep(c(0, 1), nn/2), runif(nn, 30, 70))
  beta <- c(0, 0, 5)
  epsilon <- 0

  true_lambda_probs <- happi::expit(c(xx %*% beta))
  lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)

  gamma1 <- (logit(0.8) - logit(0.2)) / (100 - 30)
  gamma0 <- logit(0.8) - gamma1 * 100

  true_f <- happi::expit(gamma0 + gamma1 * mm)

  true_prob_y_equal_1 <- true_f # correct for lambdas == 1
  true_prob_y_equal_1[lambdas == 0] <- epsilon

  ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

  hp <- happi(outcome = ys,
              covariate = xx,
              quality_var = mm,
              max_iterations = 1000,
              method = "splines",
              nstarts = 10,
              epsilon = 0,
              norm_sd = 1)
  hp %>% names
  hp$starts_df

  hp15 <- happi(outcome = ys,
                covariate = xx,
                quality_var = mm,
                max_iterations = 1000,
                method = "splines",
                nstarts = 10,
                epsilon = 0,
                norm_sd = 15)


  hp_d <- happi(outcome = ys,
                covariate = xx,
                quality_var = mm,
                method = "splines",
                epsilon = 0)

  res <- suppressWarnings(happi::happi(outcome = ys,
                                       covariate = xx,
                                       quality_var = mm,
                                       max_iterations = 1000,
                                       method = "splines",
                                       nstarts = 10,
                                       epsilon = 0,
                                       norm_sd = 1))

  # can make norm_sd larger to hit the optim error associated with infinite
  # expit(XB) for more of the starts

  res$starts_df

  # throw error if at least one start failed to converge by looking at the `starts_df`
  # data frame giving the likelihood under the null and alternative for each start
  expect_true(sum(is.na(res$starts_df)) == 0)

})
