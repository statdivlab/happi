#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("NOPE! Bad args. ARGH!", call.=FALSE)
} else {
  input_1 = as.numeric(args[1])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
  input_2 = as.numeric(args[2])
}


happi_out <- happi(outcome = ys,
                   covariate = xx,
                   quality_var = mm,
                   max_iterations = 50,
                   method="splines", firth=T, spline_df=4,
                   nstarts = 1,
                   change_threshold = 0.1,
                   epsilon = 0)
