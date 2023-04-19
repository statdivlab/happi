#!/usr/bin/env Rscript
# libraries
if (!require("argparse", quietly = TRUE))
  install.packages("argparse") # check that argparse is installed
library(argparse)

# parse data
parser <- ArgumentParser(description= 'This is a workflow for running happi')
parser$add_argument('--input', '-i', help= 'Input file')
parser$add_argument('--output', '-o', help= 'Output file')
parser$add_argument('--outcome', '-y', help= 'Outcome variable', type= 'double')
parser$add_argument('--quality_var', '-q', help= 'Quality variable', type= 'double')
parser$add_argument('--covariate', '-x', help= 'Imported variable', type= 'double')
parser$add_argument('--minit', '-a', help= 'Imported variable', type= 'double')
parser$add_argument('--maxit', '-z', help= 'Imported variable', type= 'double')
parser$add_argument('--method', '-m', help= 'Imported variable', type= 'double')
parser$add_argument('--firth', '-f', help= 'Imported variable', type= 'double')
parser$add_argument('--spline_df', '-d', help= 'Imported variable', type= 'double')
parser$add_argument('--perm', '-p', help= 'Imported variable', type= 'double')
parser$add_argument('--nstarts', '-n', help= 'Imported variable', type= 'double')
parser$add_argument('--seed', '-s', help= 'Imported variable', type= 'double')

xargs<- parser$parse_args()
# process data
my_data <- read.csv(xargs$input)
my_data <- read.csv("workflows/TM7_data.csv")

importedFactor <- xargs$myFactor
out <- dat*importedFactor

if (xargs$perm %in% c("FALSE","F")) {

  happi_out <- happi(outcome = ys,
                     covariate = xx,
                     quality_var = mm,
                     max_iterations = 50,
                     method="splines", firth=T, spline_df=4,
                     nstarts = 1,
                     change_threshold = 0.1,
                     epsilon = 0)
  
  
} else if (xargs$perm %in% c("TRUE","T")) {
  
  
}
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
