#!/usr/bin/env Rscript
# libraries
if (!require("argparse", quietly = TRUE))
  install.packages("argparse") # check that argparse is installed
library(argparse)

# parse data
parser <- ArgumentParser(description= 'This is a workflow for running happi')
parser$add_argument('--input', '-i', help= 'Input file')
parser$add_argument('--output', '-o', help= 'Output file')
parser$add_argument('--quality_var', '-q', help= 'Quality variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--myFactor', '-f', help= 'Imported variable', type= 'double')

xargs<- parser$parse_args()
# process data
dat <- readRDS(xargs$input)
importedFactor <- xargs$myFactor
out <- dat*importedFactor

# output data
write.table(out, xargs$output)

happi_out <- happi(outcome = ys,
                   covariate = xx,
                   quality_var = mm,
                   max_iterations = 50,
                   method="splines", firth=T, spline_df=4,
                   nstarts = 1,
                   change_threshold = 0.1,
                   epsilon = 0)
