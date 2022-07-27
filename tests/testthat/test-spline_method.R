library(mgcv)

test_that("method spline works", {

  data("logit_data")
  gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
  epsilon <- 0

  set.seed(506)

  replicate(10, {
    nn <- 30
    xx_sd <- 0.5
    beta <- c(1, 1)
    mm <- seq(10, 40, length.out = nn)
    true_f <- predict(gam_model, type = "response", newdata = data.frame("coverage" = mm))

    x1 <- rnorm(nn, seq(0, 1, length.out = nn), sd=xx_sd)
    xx <- cbind(1, x1)

    true_lambda_probs <- expit(c(xx %*% beta))

    lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
    true_prob_y_equal_1 <- true_f # for lambdas == 1
    true_prob_y_equal_1[lambdas == 0] <- epsilon # for lambdas == 0
    ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

    expect_type(happi(outcome = ys,
                      covariate = xx,
                      quality_var = mm,
                      max_iterations = 200,
                      nstarts = 1,
                      method = "spline",
                      change_threshold = 0.05,
                      epsilon = 0),
                "list")
  })

})
