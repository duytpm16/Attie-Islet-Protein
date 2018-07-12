### This script reads in the required data file to run the qtl2 scan1 function.
#   The proteins can be ran in 'chunks' or all at once.
#     Example for chunk size/number: 
#         There should be 5433 protein columns. If chunk size is 1000, then there should be 6 scan1 runs, 
#         with the chunk_number value being 1-6 to get the column numbers:
#         1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#         If you do not want to run in chunks, set chunk_size to 5433 and chunk_num to 1.
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
# 1: input.file:    Prefix of qtl2 input data
# 2: output.prefix: File prefix name to save the output of scan1
# 3: chunk_size:    Integer that is the chunk size.
# 4: chunk_number:  Integer that is the chunk number to run.
# 5: num_cores:     Number of cores to run
# 6: should.rankz:  Logical value to use the rankz dataset instead of normalized
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.Rdata"))
output.prefix <- args[2]
chunk_size <- as.numeric(args[3])
chunk_number <- as.numeric(args[4])
num_cores <- as.numeric(args[5])
should.rankz <- as.logical(args[6])



### Check to see if required data are loaded
stopifnot(c("norm", "pheno.dict", "genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype","display.name") %in% ls())



### If should.rankz is true used the normalized rankz dataset instead, else use the normalized dataset.
if(should.rankz){
  
  expr <- rankz
  
} else{
  
  expr <- norm
  
}




### Calculate the phenotype range to run. 
max_col = ncol(expr)
expr.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)

if(expr.rng[length(expr.rng)] > max_col) {
  
 expr.rng <- expr.rng[1]:max_col
  
} 

print(paste("Mapping columns:", expr.rng[1], "to", expr.rng[length(expr.rng)]))




### Qtl2 scan1 run
qtl <- scan1(genoprobs = genoprobs, pheno = expr[,expr.rng, drop = FALSE], kinship = K, 
            addcovar = covar, cores = num_cores)



### Save output of scan1 to current directory
saveRDS(qtl, file = paste0(output.prefix, "_chunk_", chunk_number, "_QTL.rds"))
