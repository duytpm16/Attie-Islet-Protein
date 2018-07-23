### This script reads in the required data file to run the qtl2 scan1 function.
#
### Input: 
#       1.) Path + prefix to Rdata file generated from: gather_qtl2_scan1_input.R
#       2.) Number of cores to run
#       3.) TRUE or FALSE to use rank Z data 
#
### Output: 
#       1.) Matrix of the LOD scoress for each of the phenotype that was inputted to scan1.
#
### Author: Duy Pham and Dan Gatti
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
# 1: input.file:    Prefix of qtl2 input data
# 3: num_cores:     Number of cores to run
# 4: should.rankz:  Logical value to use the rankz dataset instead of normalized
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.Rdata"))
num_cores <- as.numeric(args[2])
should.rankz <- as.logical(args[3])



### Check to see if required data are loaded
stopifnot(c("norm", "pheno.dict", "genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype","display.name") %in% ls())



### If should.rankz is true used the normalized rankz dataset instead, else use the normalized dataset.
if(should.rankz){
  
  expr <- rankz
  
} else{
  
  expr <- norm
  
}



### Qtl2 scan1 run
qtl <- scan1(genoprobs = genoprobs, pheno = expr, kinship = K, 
            addcovar = covar, cores = num_cores)


class(qtl) = c("scan1", "matrix")



### Save output of scan1 to current directory
if(should.rankz){
  saveRDS(qtl, file = paste0(args[1], "_rZ_qtl_lod.rds"))
}else{
  saveRDS(qtl, file = paste0(args[1], "qtl_lod.rds"))
}
