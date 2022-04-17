library(happi)

## yay, matches!!!

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

## pauline's estimates max iteration = 500: 
# beta0: 2.87168
# beta1: -1.521753
# LL_fullmodel: -10.48573
# LL_null: -10.72084
# LRT: 0.470219
# pvalue: 0.4928864

## estimates using happi r-package and maxit = 500 (using the same initial parameters as I used): 
## beta0 = 0.1, beta1 = 0.1, f_tilde = (1, num)
# beta0: 2.87168
# beta1: -1.521753
# LL_fullmodel: -10.47948
# LL_null: -10.71932
# LRT:  0.4796723
# pvalue: 0.4885708

happi_r27 <- happi(outcome = r27$`Ribosomal protein L27`,
                   covariate = cbind(1, r27$tongue),
                   quality_var = r27$mean_coverage,
                   max_iterations = 500,
                   nstarts = 1,
                   epsilon = 0)
### output when restored to happi's initial parameters (beta0=4, beta1 = -2, logit(f) = 0.73)
# p = 0.4671942
# LRT = 0.529
happi_r27$beta %>% tail
happi_r27$loglik %>% tail(1) 
happi_r27$f %>% tail(2) %>% logit ## doesn't match Pauline's
happi_r27$f %>% head(2) %>% logit
happi_r27$p %>% head()

happi_r27$loglik$pvalue

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
# $beta0
# [1] 4.007333
#
# $beta1
# [1] -2.985682
#
# $LRT
# [1] 5.767931
#
# $iterations_fullmodel
# [1] 100
#
# $iterations_restrictedmodel
# [1] 100
#
# $ftilde_change_fullmodel
# [1] 0.6127409
#
# $ftilde_change_restrictedmodel
# [1] 1.257598
#
# $pvalue
# [1] 0.01632124

# [1] "my_results$V3"
# [1] -8.177821
# [1] "my_restricted_results$V3"
# [1] -11.06179
### ^^ just FYI that the ftilde_change outputted pauline's code is: 
#max_{M_i} (0.001, | (f-tilde(M_i, iteration = 500) - f-tilde(M_i, iteration = 501)/(f-tilde(M_i, iteration = 500)|  )
   

names(mem)[3] <- "menbraneprotein"
happi_r27 <- happi(outcome = mem$menbraneprotein,
                   covariate = cbind(1, mem$tongue),
                   quality_var = mem$mean_coverage,
                   max_iterations = 100,
                   nstarts = 1,
                   epsilon = 0)
happi_r27$beta %>% tail # match
happi_r27$f %>% tail(1) %>% logit # not quite
happi_r27$loglik %>% tail # not quite
happi_r27$loglik$loglik %>% tail(20) %>% plot
happi_r27$loglik$loglik_null %>% tail(20) %>% plot
happi_r27$loglik$pvalue %>% tail(20) %>% plot
happi_r27$p[1, ] %>% plot
happi_r27$p[29, ] %>% plot
## yes, correct, excellent

## f's don't quite match
