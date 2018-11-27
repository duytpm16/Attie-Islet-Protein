### This script reads in the required input data file generated from gather_qtl2_scan1_input_data.R
#       to run the qtl2 scan1 function.
#   The QTL scan can be ran in 'chunks' or all at once.
#       Example for chunk size/number: 
#           Suppose there are 5433 phenotype columns. If chunk size is 1000, then there should be 6 different chunks
#             of scan1 runs, with the chunk_number value being 1-6 to get the column numbers:
#             1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#           If you do not want to run in chunks, set use_chunks to FALSE.
#
### Input:
#       1: input.file:    Path + prefix to the qtl2 input data generated from gather_qtl2_scan1_input_data.R
#       2: num_cores:     Number of cores to run
#       3: should_rankz:  Logical value to use the rankz dataset instead of normalized
#       4: use_chunks:    Logical value to run QTL scans in chunks
#       5: use_int:       Logical value to use an interaction term
#       6: chunk_number:  Numeric value of the chunk number. Not needed if use_chunks is FALSE
#       7: chunk_size:    Numeric value of chunk size. Should be consistent. Not needed if use_chunks is FALSE
#       8: int_name:      Name of the interaction term. Not needed if use_int is FALSE
#
### Output: 
#       1: Matrix containing LOD scoress for each of the phenotype that was given to scan1 at each marker.
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
options(stringsAsFactors = FALSE)
library(qtl2) 








### Command line arguments / Variables to change
# 1: input.file:    Path + prefix to the qtl2 input data generated from gather_qtl2_scan1_input_data.R
# 2: num_cores:     Number of cores to run
# 3: should_rankz:  Logical value to use the rankz dataset instead of normalized
# 4: use_chunks:    Logical value to run QTL scans in chunks
# 5: use_int:       Logical value to use an interaction term
# 6: chunk_number:  Numeric value of the chunk number. Not needed if use_chunks is FALSE
# 7: chunk_size:    Numeric value of chunk size. Should be consistent. Not needed if use_chunks is FALSE
# 8: int_name:      Name of the interaction term. Not needed if use_int is FALSE
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.RData"))
num_cores <- as.numeric(args[2])
should_rankz <- as.logical(args[3])
use_chunks <- as.logical(args[4])
use_int <- as.logical(args[5])









### Check to see if required data are loaded
stopifnot(c("norm","pheno.dict","display.name","genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype") %in% ls())










### If should_rankz is true used the normalized rankz dataset, else use the normalized dataset.
if(should_rankz){
   expr <- rankz
            
}else{
   expr <- norm
}








### If running qtl scans by chunks, get chunk size and chunk number
if(use_chunks){
   chunk_number <- as.numeric(args[6])
   chunk_size <- as.numeric(args[7])
}










### If sex should be used as an interaction term
if(use_int == TRUE & use_chunks == TRUE){
   int_name <- args[8]
   formula <- as.formula(paste0('~', int_name))
   int_mat <- model.matrix(formula, data = samples)[,-1]
}

if(use_int == TRUE & use_chunks == FALSE){
   int_name <- args[6]
   formula <- as.formula(paste0('~', int_name))
   int_mat <- model.matrix(formula, data = samples)[,-1]                 
}

if(use_int == FALSE){
   int_mat <- NULL
}









### Running QTL2 scan1 function
if(use_chunks == FALSE){
            
   # qtl2 scan1 without chunks
   qtl <- scan1(genoprobs = genoprobs, 
                pheno = expr,
                kinship = K, 
                addcovar = covar,
                intcovar = int_mat, 
                cores = num_cores)
    
   # Set class of qtl matrix
   class(qtl) = c("scan1", "matrix")

            
   # Save output of scan1 to current directory
   if(use_int){
      if(should_rankz){
         saveRDS(qtl, file = paste0(args[1], "_",int_name,"_int_rZ_qtl_lod.rds"))
      }else{
         saveRDS(qtl, file = paste0(args[1], "_",int_name,"_int_norm_qtl_lod.rds"))
      }
   }else{
      if(should_rankz){
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

   }
   print(paste("Mapping columns:", pheno.rng[1], "to", pheno.rng[length(pheno.rng)]))

            
            
   # Run QTL scan with chunks
   qtl <- scan1(genoprobs = genoprobs,
                pheno = expr[,pheno.rng,drop = FALSE],
                kinship = K,
                addcovar = covar,
                intcovar = int_mat,
                cores = num_cores)

            
            
            
            
   # Save the QTL chunk scan
   if(use_int){
      if(should_rankz){
         saveRDS(qtl, file = paste0(args[1], "_",int_name,"_int_rZ_qtl_lod_chunk_",chunk_numer,".rds"))
      
      }else{
         saveRDS(qtl, file = paste0(args[1], "_",int_name,"_int_norm_qtl_lod_chunk_",chunk_number,".rds"))
      }
     
   }else{
      if(should_rankz){
         saveRDS(qtl, file = paste0(args[1], "_rZ_qtl_lod_chunk_",chunk_number,".rds"))
      
      }else{
         saveRDS(qtl, file = paste0(args[1], "_norm_qtl_lod_chunk_",chunk_number,".rds"))
      }
   }		   
}
