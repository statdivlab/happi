library(mgcv)
test_that("likelihood is higher under alternative than null", {

  data("logit_data")
  gam_model <- mgcv::gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
  epsilon <- 0

  set.seed(3) # or 5
  # find_failures <- function(seed) {
  #   set.seed(seed)
  nn <- 100
  xx_sd <- 0.5
  beta <- c(0, 0)
  mm <- seq(10, 40, length.out = nn)
  true_f <- stats::predict(gam_model, type = "response", newdata = data.frame("coverage" = mm))

  x1 <- rnorm(nn, seq(0, 1, length.out = nn), sd=xx_sd)
  xx <- cbind(1, x1)

  true_lambda_probs <- happi::expit(c(xx %*% beta))

  lambdas <- rbinom(n = nn, size = 1, prob = true_lambda_probs)
  true_prob_y_equal_1 <- true_f # for lambdas == 1
  true_prob_y_equal_1[lambdas == 0] <- epsilon # for lambdas == 0
  ys <- rbinom(n = nn, size = 1, prob = true_prob_y_equal_1)

  happi_spline_firth <- happi::happi(outcome = ys,
                              covariate = xx,
                              quality_var = mm,
                              max_iterations = 1000,
                              nstarts = 1,
                              method = "splines",
                              firth=T,
                              change_threshold = 0.1,
                              epsilon = 0)
  happi_spline_logreg <- happi::happi(outcome = ys,
                               covariate = xx,
                               quality_var = mm,
                               max_iterations = 1000,
                               nstarts = 1,
                               method = "splines",
                               firth=F,
                               change_threshold = 0.1,
                               epsilon = 0)

   happi_isotone_firth <- happi::happi(outcome = ys,
                               covariate = xx,
                               quality_var = mm,
                               max_iterations = 1000,
                               nstarts = 1,
                               method = "isotone",
                               firth=T,
                               change_threshold = 0.1,
                               epsilon = 0)
   happi_isotone_logreg <- happi::happi(outcome = ys,
                              covariate = xx,
                                quality_var = mm,
                                max_iterations = 1000,
                                nstarts = 1,
                                method = "isotone",
                                firth=F,
                                change_threshold = 0.1,
                                epsilon = 0)


  expect_gte(tail(happi_spline_firth$loglik$LRT_nopenalty[!is.na(happi_spline_firth$loglik$LRT_nopenalty)], 1), 0)
  expect_gte(tail(happi_spline_logreg$loglik$LRT_nopenalty[!is.na(happi_spline_logreg$loglik$LRT_nopenalty)], 1), 0)
  expect_gte(tail(happi_isotone_firth$loglik$LRT_nopenalty[!is.na(happi_isotone_firth$loglik$LRT_nopenalty)], 1), 0)
  expect_gte(tail(happi_isotone_logreg$loglik$LRT_nopenalty[!is.na(happi_isotone_logreg$loglik$LRT_nopenalty)], 1), 0)

})
