###############################################################################################################################
#
#
#   This script finds LOD peaks above a given threshold 
#      and computes the allele effects at each QTL for additive scans only.
#
#
#
#   Input:
#       1.) Path + prefix to .RData file generarted from: qtl2_gather_scan1_input.R
#       2.) Numeric LOD threshold used to find significant QTLs.
#       3.) Number of cores to run
#       4.) TRUE or FALSE value to use rank Z data.
#       5.) TRUE or FALSE to use interaction lod scan data.
#       6.) Name of the interaction term used. Not needed if 5 is FALSE
#
#
#
#   Output:
#       1.) .rds file containing the lod peaks above threshold and allele effects at QTLS (Only if 5 in Input is FALSE).
#
#
#   Author: Duy Pham
#   Date:   July 11, 2018
#   E-mail: duy.pham@jax.org
#
###############################################################################################################################

# Load libraries
library(qtl2, lib.loc = '~/Rlibs')

# Setting options
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Variables to change
#       1.) Path + prefix to .RData file generarted from: qtl2_gather_scan1_input.R
#       2.) Numeric LOD threshold used to find significant QTLs.
#       3.) Number of cores to run
#       4.) TRUE or FALSE value to use rank Z data.
#       5.) TRUE or FALSE to use interaction lod scan data.
#       6.) Name of the interaction term used. Not needed if 5 is FALSE
args = commandArgs(trailingOnly = TRUE)

prefix <- args[1]
threshold <- as.numeric(args[2])
num_cores <- as.numeric(args[3])
should_rankz <- as.logical(args[4]) 
use_int <- as.logical(args[5])








### Load in the data
if(use_int){
   # If using interaction scan
   int_name <- args[6]

   load(paste0(prefix,"_qtl2_input.RData"))
  
   stopifnot(c("norm","pheno.dict","display.name","genoprobs", "K", "map", "markers", "covar", 
               "covar.factors", "samples", "rankz", "raw","datatype") %in% ls())
  
   if(should_rankz){
      scan1_data <- readRDS(paste0(prefix,"_",int_name,"_int_rZ_qtl_lod.rds"))
      expr <- rankz
   }else{
      scan1_data <- readRDS(paste0(prefix,"_",int_name,"_int_norm_qtl_lod.rds"))
      expr <- norm
   }
  
  
}else{
   # If using additive scan
   load(paste0(prefix,"_qtl2_input.Rdata"))
    
   # Check to see if required QTL2 data are loaded
   stopifnot(c("norm","pheno.dict","display.name","genoprobs", "K", "map", "markers", "covar", 
               "covar.factors", "samples", "rankz", "raw","datatype") %in% ls())
  
  
  if(should_rankz){
     scan1_data <- readRDS(paste0(prefix,"_rZ_qtl_lod.rds"))
     expr <- rankz
  }else{
     scan1_data <- readRDS(paste0(prefix,"_norm_qtl_lod.rds"))
     expr <- norm
  }
}









### Run the find_peaks function and stores the output to a variable called "lod.peaks"
lod.peaks <- find_peaks(scan1_output =  scan1_data, map = map, threshold = threshold, cores = num_cores)








### Setting up the lod.peaks dataframe in QTL Viewer format.
marker.id <- paste0(as.character(lod.peaks$chr), '_', round(lod.peaks$pos * 1000000))
annot.id <- lod.peaks[,'lodcolumn']
lod.peaks <- cbind(annot.id, marker.id, lod.peaks[,c('lod','chr','pos')])
colnames(lod.peaks) <- c('annot.id','marker.id','lod','qtl.chr','qtl.pos')
lod.peaks$qtl.chr <- as.character(lod.peaks$qtl.chr)










### BLUP scan on additive scan
if(use_int == FALSE){
  
  
   # Adding 8 more columns to lod.peaks for BLUP mapping on additive scans only
   lod.peaks = cbind(lod.peaks, matrix(0, nrow = nrow(lod.peaks), ncol = 8, 
                                       dimnames = list(NULL, LETTERS[1:8])))
  
  
  
   # BLUP mapping begin
   for(i in 1:nrow(lod.peaks)) {
  
       chr  = lod.peaks$qtl.chr[i]
       mkr  = lod.peaks$marker.id[i]
       gene = lod.peaks$annot.id[i]
  
   # BLUP scan at QTL.
       gp = genoprobs[,chr]
       gp[[1]] = gp[[1]][,,mkr,drop=FALSE]
       blup = scan1blup(genoprobs = gp, pheno = expr[,gene, drop = FALSE],
                        kinship = K[[chr]], addcovar = covar, cores = num_cores)
      
       lod.peaks[i,6:13] = blup[1,1:8]
     
   } # for(i)
}













# Save the output of the find_peaks function
if(use_int){
  
  if(should_rankz){
    saveRDS(lod.peaks,file = paste0(prefix,"_",int_name,"_int_rZ_lodpeaks_",threshold,".rds"))
  else{
    saveRDS(lod.peaks,file = paste0(prefix,"_",int_name,"_int_norm_lodpeaks_",threshold,".rds"))
  }
    
}else{
    
  if(should_rankz){
    saveRDS(lod.peaks,file = paste0(prefix,"_rZ_lodpeaks_",threshold,".rds"))
  }else{
    saveRDS(lod.peaks,file = paste0(prefix,"_norm_lodpeaks_",threshold,".rds"))
  }
    
}
