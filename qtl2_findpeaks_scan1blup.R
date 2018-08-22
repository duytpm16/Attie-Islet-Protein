### This script finds LOD peaks above a given threshold 
#     and computes the allele effects.
#
### Input:
#       1.) Path + prefix to .RData file generarted from: qtl2_gather_scan1_input.R and .rds file from qtl2_scan1.R.
#             Asssuming both are in the same directory.
#       2.) Numeric threshold used to find QTLs with LODs score above.
#       3.) Number of cores to run
#       4.) TRUE or FALSE value to use rank Z data.
#
### Output:
#       1.) RDS data containing the lod peaks above threshold and allele effects.
#
### Author: Duy Pham
### Date:   July 11, 2018
### E-mail: duy.pham@jax.org
################################################################################

# Load libraries
library(qtl2, lib.loc = '~/Rlibs')

# Setting options
options(stringsAsFactors = FALSE)
options(scipen = 999)



### Variables to change
#     1. Prefix to .RData file containing: 
#             norm, raw, rankz, covar, covar.factors, genoprobs, samples, map, marker_map, K, datatype,  and pheno.dict
#         AND 
#             large scan1 LOD matrix. Assuming they are stored in the same directory
#     2. Threshold value for LOD peaks.
#     3. Number of cores to run
#     4. TRUE or FALSE value to use the rankz data
#     5. TRUE or FALSE value to use the scan1 matrix with sex as an interaction term
args = commandArgs(trailingOnly = TRUE)

prefix <- args[1]
threshold <- as.numeric(args[2])
num_cores <- as.numeric(args[3])
should_rankz <- as.logical(args[4])
use_sexint <- as.logical(args[5])



### Load in the data
if(should_rankz){
  
  load(paste0(prefix,"_qtl2_input.Rdata"))
  
  if(use_sexint){
    scan1_data <- readRDS(paste0(prefix,"_sexint_rZ_qtl_lod.rds"))
  }else{
    scan1_data <- readRDS(paste0(prefix,"_rZ_qtl_lod.rds"))
  }
  
}else{
  
  load(paste0(prefix,"_qtl2_input.Rdata"))
  
  if(use_sexint){
    scan1_data <- readRDS(paste0(prefix,"_sexint_norm_qtl_lod.rds"))
  }else{
    scan1_data <- readRDS(paste0(prefix,"_norm_qtl_lod.rds"))
  }
}



### Check to see if required QTL2 data are loaded
stopifnot(c("norm", "pheno.dict", "genoprobs", "K", "map", "markers", "covar", 
            "covar.factors", "samples", "rankz", "raw","datatype","display.name") %in% ls()
)

if(should_rankz){
     expr <- rankz
}else{
     expr <- norm
}	


### Run the find_peaks function and stores the output to a variable called "lod.peaks"
lod.peaks <- find_peaks(scan1_output =  scan1_data, map = map, threshold = threshold, cores = num_cores)



### Setting up the lod.peaks dataframe in QTL Viewer format.
marker.id <- paste0(as.character(lod.peaks$chr), '_', round(lod.peaks$pos * 1000000))

annot.id <- lod.peaks[,'lodcolumn']
lod.peaks <- cbind(annot.id, marker.id, lod.peaks[,c('lod','chr','pos')])
lod.peaks$chr <- as.character(lod.peaks$chr)



### Adding 8 more columns to lod.peaks for BLUP mapping
lod.peaks = cbind(lod.peaks, matrix(0, nrow = nrow(lod.peaks), ncol = 8, 
                                    dimnames = list(NULL, LETTERS[1:8])))



# BLUP mapping.
if(use_sexint == FALSE){
    for(i in 1:nrow(lod.peaks)){

      chr  = lod.peaks$chr[i]
      mkr  = lod.peaks$marker.id[i]
      gene = lod.peaks$annot.id[i]

      # Scans.
      gp = genoprobs[,chr]
      gp[[1]] = gp[[1]][,,mkr,drop=FALSE]
      blup = scan1blup(genoprobs = gp, pheno = expr[,gene, drop = FALSE],
                       kinship = K[[chr]], addcovar = covar, cores = num_cores)
      lod.peaks[i,6:13] = blup[1,1:8]
    }
}


# Save the output of the find_peaks function
rm(scan1_data, annot.id, marker.id, i, mkr, gene, chr, gp, blup)

if(should_rankz){
  
  if(use_sexint){
    saveRDS(lod.peaks,file = paste0(prefix,'_sexint_rZ_lod_peaks_',threshold,'.rds'))
  else{
    saveRDS(lod.peaks,file = paste0(prefix,'_rZ_lod_peaks_',threshold,'.rds'))
  }
    
}else{
    
  if(use_sexint){
    saveRDS(lod.peaks,file = paste0(prefix,'_sexint_norm_lod_peaks_',threshold,'.rds'))
  }else{
    saveRDS(lod.peaks,file = paste0(prefix,'_norm_lod_peaks_',threshold,'.rds'))
  }
    
}
