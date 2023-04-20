#!/usr/bin/env Rscript
# libraries
if (!require("argparse", quietly = TRUE))
  install.packages("argparse") # check that argparse is installed
library(argparse)
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr") # check that dplyr is installed
library(dplyr)
if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr") # check that magrittr is installed
library(magrittr)
if (!require("parallel", quietly = TRUE))
  install.packages("parallel") # check that parallel is installed
library(parallel)
if (!require("happi", quietly = TRUE))
  install.packages("happi") # check that parallel is installed
library(happi)

# parse data
parser <- ArgumentParser(description= 'This is a workflow for running happi')
parser$add_argument('--input', '-i', help= 'Input file containing gene presence absence')
parser$add_argument('--metadata', '-m', help= 'metadata file')
parser$add_argument('--output', '-o', help= 'output file')
parser$add_argument('--quality_var', '-q', help= 'quality variable')
parser$add_argument('--covariate', '-x', help= 'covariate variable')
parser$add_argument('--epsilon', '-e', help= 'epsilon variable', type= 'integer')
parser$add_argument('--minit', '-a', help= 'minimum iterations', type= 'integer')
parser$add_argument('--maxit', '-z', help= 'maximum iterations', type= 'integer')
parser$add_argument('--method', '-M', help= 'splines or isotone')
parser$add_argument('--firth', '-f', help= 'firth penalty')
parser$add_argument('--spline_df', '-df', help= 'splines degrees of freedom', type= 'integer')
parser$add_argument('--perm', '-p', help= 'use permutation test')
parser$add_argument('--num_perm', '-B', help= 'number of permutations', type= 'integer')
parser$add_argument('--nstarts', '-n', help= 'number of initial starts', type= 'integer')
parser$add_argument('--seed', '-s', help= 'seed number', type= 'integer')
parser$add_argument('--change_threshold', '-c', help= 'change threshold for E-M', type= 'double')
parser$add_argument('--cores', '-t', help= 'number of cores', type= 'integer')

 xargs <- parser$parse_args()
 print(xargs)
# process data
 my_gene_presence <- read.csv(xargs$input)
 my_metadata <- read.csv(xargs$metadata)
 
 my_gene_presence <- my_gene_presence %>% dplyr::arrange(ID)
 my_metadata <- my_metadata %>% dplyr::arrange(ID)
 # Ordering everything to make sure input data aligns 
 if (dim(my_gene_presence)[1] != dim(my_metadata)[1]){
   message("The number of samples in your gene presence file do not match your metadata file. Please check the format of your data.")
 }
 
 my_quality_var <- my_metadata %>% as.data.frame() %>% select(as.name(xargs$quality_var))
 #my_quality_var <- my_metadata %>% as.data.frame() %>% select(as.name("mean_coverage"))
 
 covariate <- as.name(xargs$covariate)
 
 x_matrix <- stats::model.matrix(~eval(covariate), data = my_metadata)
 
 my_gene_parallel <- my_gene_presence %>% dplyr::select(-ID)
 
 firth <- xargs$firth
 perm <- xargs$perm 
 
 run_happi_parallel_perm_true_firth_true <- function(colnum) {
   happi_results <- happi::happi(outcome=unlist(my_gene_parallel[,colnum]), 
                                 covariate = x_matrix, 
                                 quality_var =  my_metadata[,xargs$quality_var],
                                 method = xargs$method,
                                 firth = T, 
                                 spline_df = as.numeric(xargs$spline_df),
                                 max_iterations = as.numeric(xargs$maxit), 
                                 min_iterations = as.numeric(xargs$minit), 
                                 change_threshold = as.numeric(xargs$change_threshold), 
                                 epsilon = as.numeric(xargs$epsilon), 
                                 seed = as.numeric(xargs$seed), 
                                 nstarts = as.numeric(xargs$nstarts))
   npLRT_results <- happi::npLRT(happi_results, 
                          firth = T, 
                          spline_df = as.numeric(xargs$spline_df),
                          max_iterations = as.numeric(xargs$maxit), 
                          min_iterations = as.numeric(xargs$minit), 
                          change_threshold = as.numeric(xargs$change_threshold), 
                          epsilon = as.numeric(xargs$epsilon),
                          nstarts = as.numeric(xargs$nstarts), 
                          P = as.numeric(xargs$num_perm))
   
   return(npLRT_results)
 }
 run_happi_parallel_perm_true_firth_false <- function(colnum) {
   happi_results <- happi::happi(outcome=unlist(my_gene_parallel[,colnum]), 
                                 covariate = x_matrix, 
                                 quality_var =  my_metadata[,xargs$quality_var],
                                 method = xargs$method,
                                 firth = F, 
                                 spline_df = as.numeric(xargs$spline_df),
                                 max_iterations = as.numeric(xargs$maxit), 
                                 min_iterations = as.numeric(xargs$minit), 
                                 change_threshold = as.numeric(xargs$change_threshold), 
                                 epsilon = as.numeric(xargs$epsilon), 
                                 seed = as.numeric(xargs$seed), 
                                 nstarts = as.numeric(xargs$nstarts))
   npLRT_results <- happi::npLRT(happi_results, 
                          firth = F, 
                          spline_df = as.numeric(xargs$spline_df),
                          max_iterations = as.numeric(xargs$maxit), 
                          min_iterations = as.numeric(xargs$minit), 
                          change_threshold = as.numeric(xargs$change_threshold), 
                          epsilon = as.numeric(xargs$epsilon), 
                          nstarts = as.numeric(xargs$nstarts), 
                          P = as.numeric(xargs$num_perm))
   return(npLRT_results)
 }
 
 run_happi_parallel_firth_true <- function(colnum) {
   happi_results <- happi::happi(outcome=unlist(my_gene_parallel[,colnum]), 
                                 covariate = x_matrix, 
                                 quality_var =  my_metadata[,xargs$quality_var],
                                 method = xargs$method,
                                 firth = T, 
                                 spline_df = as.numeric(xargs$spline_df),
                                 max_iterations = as.numeric(xargs$maxit), 
                                 min_iterations = as.numeric(xargs$minit), 
                                 change_threshold = as.numeric(xargs$change_threshold), 
                                 epsilon = as.numeric(xargs$epsilon), 
                                 seed = as.numeric(xargs$seed), 
                                 nstarts = as.numeric(xargs$nstarts))
   return(happi_results)
 }
 
   run_happi_parallel_firth_false <- function(colnum) {
     happi_results <- happi::happi(outcome=unlist(my_gene_parallel[,colnum]), 
                                   covariate = x_matrix, 
                                   quality_var =  my_metadata[,xargs$quality_var],
                                   method = xargs$method,
                                   firth = F, 
                                   spline_df = as.numeric(xargs$spline_df),
                                   max_iterations = as.numeric(xargs$maxit), 
                                   min_iterations = as.numeric(xargs$minit), 
                                   change_threshold = as.numeric(xargs$change_threshold), 
                                   epsilon = as.numeric(xargs$epsilon), 
                                   seed = as.numeric(xargs$seed), 
                                   nstarts = as.numeric(xargs$nstarts))
     return(happi_results)
  }

 
 run_happi_parallel_perm_true_firth_true <- function(colnum) {
     happi_results <- happi::happi(outcome=unlist(my_gene_parallel[,colnum]), 
                                   covariate = x_matrix, 
                                   quality_var =  my_metadata[,xargs$quality_var],
                                   method = xargs$method,
                                   firth = T, 
                                   spline_df = as.numeric(xargs$spline_df),
                                   max_iterations = as.numeric(xargs$maxit), 
                                   min_iterations = as.numeric(xargs$minit), 
                                   change_threshold = as.numeric(xargs$change_threshold), 
                                   epsilon = as.numeric(xargs$epsilon), 
                                   seed = as.numeric(xargs$seed), 
                                   nstarts = as.numeric(xargs$nstarts))
     npLRT_results <- happi::npLRT(happi_results, 
                                   firth = T, 
                                   spline_df = as.numeric(xargs$spline_df),
                                   max_iterations = as.numeric(xargs$maxit), 
                                   min_iterations = as.numeric(xargs$minit), 
                                   change_threshold = as.numeric(xargs$change_threshold), 
                                   epsilon = as.numeric(xargs$epsilon),
                                   nstarts = as.numeric(xargs$nstarts), 
                                   P = as.numeric(xargs$num_perm))
     
     return(npLRT_results)
   }
 if (perm == "TRUE" & firth == "TRUE") {
   
  set.seed(as.numeric(xargs$seed))
  my_results <- parallel::mclapply(1:length(my_gene_parallel), run_happi_parallel_perm_true_firth_true, mc.cores = as.numeric(xargs$cores))
  pvalue_happi <- my_results %>%  unlist
 
  } else if (perm == "TRUE" & firth == "FALSE"){
    
   set.seed(as.numeric(xargs$seed))
   my_results <- parallel::mclapply(1:length(my_gene_parallel), run_happi_parallel_perm_true_firth_false, mc.cores = as.numeric(xargs$cores))
   pvalue_happi <- my_results %>%  unlist
   
   } else if (perm == "FALSE" & firth == "TRUE"){
   
   set.seed(as.numeric(xargs$seed))
   my_results <- parallel::mclapply(1:length(my_gene_parallel), run_happi_parallel_firth_true, mc.cores = as.numeric(xargs$cores))
   pvalue_happi <- lapply(my_results, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist

   } else if (perm == "FALSE" & firth == "FALSE"){
     
   set.seed(as.numeric(xargs$seed))
   my_results <- parallel::mclapply(1:length(my_gene_parallel), run_happi_parallel_firth_false, mc.cores = as.numeric(xargs$cores))
   pvalue_happi <- lapply(my_results, function(x) tail(x$loglik$pvalue_nopenalty[!is.na(x$loglik$pvalue_nopenalty)], 1)) %>% unlist
   
 }
   
 happi_results <- tibble("gene" = colnames(my_gene_parallel)[1:length(my_gene_parallel)], 
                          pvalue_happi) 
 
 write.csv(happi_results, xargs$output, row.names = F)
 