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



### Command line arguments.
# 1: input.file:    Prefix of qtl2 input data
# 3: num_cores:     Number of cores to run
# 4: should.rankz:  Logical value to use the rankz dataset instead of normalized
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.RData"))
num_cores <- as.numeric(args[2])
should.rankz <- as.logical(args[3])
use.diff <- as.logical(args[4])
use_chunks <- as.logical(args[5])

if(use_chunks){
  chunk_number <- as.numeric(args[6])
  chunk_size <- as.numeric(args[7])
}


### Check to see if required data are loaded
stopifnot(c("genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype") %in% ls())



### If should.rankz is true used the normalized rankz dataset instead, else use the normalized dataset.
if(should.rankz){
  
  expr <- rankz
  
} else{
  
  expr <- norm
  
}

if(use.diff){
 int_covar_name <- args[5]
 formula <- as.formula(paste0('~', int_covar_name))
 int_mat <- model.matrix(formula, data = samples)[,-1]
}else{
 int_mat <- NULL
}

if(use_chunks == FALSE){
  
    ### Qtl2 scan1 run without chunks
    qtl <- scan1(genoprobs = genoprobs, 
                 pheno = expr,
                 kinship = K, 
                 addcovar = covar,
                 intcovar = int_mat, 
                 cores = num_cores)
    
    class(qtl) = c("scan1", "matrix")

    ### Save output of scan1 to current directory
    if(use.diff){
       if(should.rankz){
          saveRDS(qtl, file = paste0(args[1], "_sexint_rZ_qtl_lod.rds"))
       }else{
          saveRDS(qtl, file = paste0(args[1], "_sexint_qtl_lod.rds"))
       }
    }else{
       if(should.rankz){
          saveRDS(qtl, file = paste0(args[1], "_rZ_qtl_lod.rds"))
       }else{
          saveRDS(qtl, file = paste0(args[1], "_qtl_lod.rds"))
       }
    }
    
}else{
   
   # Get chunk range
   max_col = ncol(expr)
   pheno.rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
   if(pheno.rng[length(pheno.rng)] > max_col) {

      pheno.rng = pheno.rng[1]:max_col

   } # if(pheno.rng[length(pheno.rng)] > max_col)

   print(paste("Mapping columns:", pheno.rng[1], "to", pheno.rng[length(pheno.rng)]))

   # Run QTL scan with chunk
   qtl <- scan1(genoprobs = genoprobs,
             pheno = expr[,pheno.rng,drop = FALSE],
             kinship = K,
             addcovar = covar,
             intcovar = int_mat,
             cores = num_cores)

   # Save the QTL chunk scan
   if(use.diff){
      if(should.rankz){
         saveRDS(qtl, file = paste0(args[1], "_sexint_rZ_qtl_lod_",chunk_numer,".rds"))
      }else{
         saveRDS(qtl, file = paste0(args[1], "_sexint_qtl_lod_",chunk_number,".rds"))
      }
   }else{
      if(should.rankz){
         saveRDS(qtl, file = paste0(args[1], "_rZ_qtl_lod_",chunk_number,".rds"))
      }else{
         saveRDS(qtl, file = paste0(args[1], "_qtl_lod_",chunk_number,".rds"))
      }
   }		   
}