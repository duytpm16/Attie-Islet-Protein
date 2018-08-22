### This script reads in the data file generated from gather_qtl2_scan1_input_data.R
#      to run the qtl2 scan1 function.
#
### Input: 
#      1.) Path + prefix to .RData file generated from gather_qtl2_scan1_input_data.R
#      2.) Number of cores to run
#      3.) TRUE or FALSE to use the rankz expression matrix
#      4.) TRUE or FALSE to use sex as an interaction term
#
### Output: 
#      1.) A scan1 matrix of LOD scores for each of the phenotype/traits at each marker
#
### Author: Duy Pham
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org
####################################################################################################################



### Install required library packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl", "dplyr))
# library(devtools)
# install_github("rqtl/qtl2", "rqtl/qtl2convert")



### Load required library packages
library(qtl2, lib.loc = '~/Rlibs') 
library(qtl2convert,lib.loc = '~/Rlibs')
library(RSQLite, lib.loc = '~/Rlibs')



### Command line arguments.
# 1: input.file:     Prefix of qtl2 input data
# 2: num_cores:      Number of cores to run
# 3: should_rankz:   Logical value to use the rankz dataset instead of normalized
# 4: use_sexint:     Logical value to use sex as an interaction term
# 5: int_covar_name: Name of the column for sex in samples dataframe (Line 70). Don't need to specify if 4 is FALSE
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.Rdata"))
num_cores <- as.numeric(args[2])
should_rankz <- as.logical(args[3])
use_sexint <- as.logical(args[4])



### Check to see if required data are loaded
stopifnot(c("norm", "pheno.dict", "genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype","display.name") %in% ls())



### If should.rankz is true used the normalized rankz dataset instead, else use the normalized dataset.
if(should_rankz){
  
  expr <- rankz
  
} else{
  
  expr <- norm
  
}



### If use.diff (sex interaction) create a new model.matrix for sex which will be used for the intcovar parameter in scan1
if(use_diff){
 int_covar_name <- args[5]
 formula <- as.formula(paste0('~', int_covar_name))
 int_mat <- model.matrix(formula, data = samples)[,-1]
}else{
 int_mat <- NULL
}



### Qtl2 scan1 run
qtl <- scan1(genoprobs = genoprobs, 
             pheno = expr,
             kinship = K, 
             addcovar = covar,
             intcovar = int_mat, 
             cores = num_cores)


class(qtl) = c("scan1", "matrix")


### Save output of scan1 to current directory
if(use_diff){
   if(should_rankz){
      saveRDS(qtl, file = paste0(args[1], "_sexint_rZ_qtl_lod.rds"))
   }else{
      saveRDS(qtl, file = paste0(args[1], "_sexint_qtl_lod.rds"))
   }
}else{
   if(should_rankz){
      saveRDS(qtl, file = paste0(args[1], "_rZ_qtl_lod.rds"))
   }else{
      saveRDS(qtl, file = paste0(args[1], "_qtl_lod.rds"))
   }
}   
