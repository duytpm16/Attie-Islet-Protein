### This script reads in the required data file to run the qtl2 scan1 function.
#
#
### Input: qtl2 data file generated from:
#
### Output: Dataframe of the LOD scoress for each of the proteins that was given to scan1.
#
### Author: Duy Pham, 'phenotype range run' was taken from Dan Gatti
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org
####################################################################################################################



### Install required library packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl", "dplyr))
# library(devtools)
# install_github("rqtl/qtl2", "rqtl/qtl2convert")



### Load required library packages
library(qtl2) 
library(qtl2convert)
library(RSQLite)
library(dplyr)



### Command line arguments.
# 1: input.file:    Prefix of qtl2 input data. If the file is not in the same directory use path to dir and prefix.
# 2: output.prefix: File prefix name to save the output of scan1
# 3: num_cores:     Number of cores to run
# 4: should.rankz:  Logical value to use the rankz dataset instead of normalized
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.Rdata"))
output.prefix <- args[2]
num_cores <- as.numeric(args[3])
should.rankz <- as.logical(args[4])



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
qtl <- scan1(genoprobs = genoprobs, pheno = expr[,expr.rng, drop = FALSE], kinship = K, 
            addcovar = covar, cores = num_cores)



### Save output of scan1 to current directory
saveRDS(qtl, file = paste0(output.prefix, "_chunk_", chunk_number, "_QTL.rds"))
