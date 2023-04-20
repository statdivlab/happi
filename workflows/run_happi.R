#!/usr/bin/env Rscript
# libraries
if (!require("argparse", quietly = TRUE))
  install.packages("argparse") # check that argparse is installed
library(argparse)
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr") # check that argparse is installed
library(dplyr)
if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr") # check that argparse is installed
library(magrittr)
if (!require("parallel", quietly = TRUE))
  install.packages("parallel") # check that argparse is installed
library(parallel)


# parse data
parser <- ArgumentParser(description= 'This is a workflow for running happi')
parser$add_argument('--input', '-i', help= 'Input file containing gene presence absence')
parser$add_argument('--metadata', '-t', help= 'metadata file')
parser$add_argument('--output', '-o', help= 'output file')
parser$add_argument('--quality_var', '-q', help= 'quality variable', type= 'double')
parser$add_argument('--covariate', '-x', help= 'covariate variable', type= 'double')
parser$add_argument('--epsilon', '-e', help= 'epsilon variable', type= 'integer')
parser$add_argument('--minit', '-a', help= 'minimum iterations', type= 'integer')
parser$add_argument('--maxit', '-z', help= 'maximum iterations', type= 'integer')
parser$add_argument('--method', '-m', help= 'splines or isotone', type= 'double')
parser$add_argument('--firth', '-f', help= 'firth penalty', type= 'double')
parser$add_argument('--spline_df', '-d', help= 'splines degrees of freedom', type= 'integer')
parser$add_argument('--perm', '-p', help= 'use permutation test', type= 'double')
parser$add_argument('--num_perm', '-B', help= 'number of permutations', type= 'integer')
parser$add_argument('--nstarts', '-n', help= 'number of initial starts', type= 'integer')
parser$add_argument('--seed', '-s', help= 'seed number', type= 'integer')
parser$add_argument('--change_threshold', '-c', help= 'change threshold for E-M', type= 'integer')
parser$add_argument('--cores', '-q', help= 'number of cores', type= 'integer')

xargs<- parser$parse_args()
# process data
my_gene_presence <- read.csv(xargs$input)
my_metadata <- read.csv(xargs$metadata) 
my_gene_presence <- tibble::tibble(read.csv("workflows/TM7_genes_presence_table.csv"))
my_metadata <- tibble::tibble(read.csv("workflows/TM7_metadata.csv"))

my_gene_presence <- my_gene_presence %>% dplyr::arrange(ID)
my_metadata <-my_metadata %>% dplyr::arrange(ID)

if (dim(my_gene_presence)[1] != dim(my_metadata)[1]){
  message("The number of samples in your gene presence file do not match your metadata file. Please check the format of your data.")
}
quality_var <- xargs$quality_var

x_matrix <- model.matrix(~xargs$covariate, data = my_metadata)

my_gene_parallel <- my_gene_presence %>% dplyr::select(-ID)

run_happi_parallel <- function(colnum) {
  happi_results <- happi(outcome=unlist(my_gene_presence[,colnum]), 
                         covariate = x_matrix, 
                         quality_var = my_metadata$quality_var,
                         method = xargs$method, 
                         firth = xargs$firth, 
                         spline_df = as.numeric(xargs$spline_df),
                         max_iterations = as.numeric(xargs$maxit), 
                         min_iterations = as.numeric(xargs$minit), 
                         change_threshold = as.numeric(xargs$change_threshold), 
                         epsilon = as.numeric(xargs$epsilon), 
                         seed = as.numeric(xargs$seed), 
                         nstarts = as.numeric(xargs$nstarts))
}
set.seed(as.numeric(xargs$seed))
my_results <- mclapply(1:length(my_gene_parallel), run_happi_parallel, mc.cores = as.numeric(xargs$cores))
pvalue_happi_a <- lapply(my_results, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist

happi_results <- tibble("gene" = colnames(my_gene_parallel)[1:length(my_gene_parallel)], 
                         pvalue_happi_a) 

write.table(happi_results, xargs$output)

