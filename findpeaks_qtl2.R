### Returns a find_peaks dataframe consisting of LOD peaks above a given threshold and allele effects.
#
### Input: 
#      1.) Prefix to qtl2 input RData. If not in the directory, use full path to directory and prefix
#      2.) Genoprobs file
#      3.) Threshold value
#      4.) Number of cores to run on
#      5.) TRUE or FALSE value to use the rankz data
#
### Output:
#      1.) rds file that contains a dataframe with lod peaks above the given threshold and allele effects.
#
### Author: Duy Pham and Dan Gatti
### Date:   July 11, 2018
### E-mail: duy.pham@jax.org
################################################################################

# Load libraries
library(qtl2)

# Setting options
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Variables to change
#     1. Prefix to Rdata file containing: 
#             norm, raw, rankz, covar, covar.factors, genoprobs, samples, map, marker_map, K, datatype,  and pheno.dict
#         AND 
#             large scan1 LOD matrix. Assuming they are store in the same directory
#     2. Threshold value for LOD peaks.
#args = commandArgs(trailingOnly = TRUE)

prefix = args[1]
genoprobs <- readRDS(args[2])
threshold = readRDSic(args[3])
num_cores = as.numeric(args[4])
should_rankz = as.logical(args[5])


### Load in the data
qtl2_data <- load(paste0(prefix,"_qtl2_input.Rdata"))
scan1_data <- readRDS(paste0(prefix,"_qtl_lod.rds"))



### Check to see if required QTL2 data are loaded
stopifnot(c("norm", "pheno.dict", "genoprobs", "K", "map", "markers", "covar",
           "covar.factors", "samples", "rankz", "raw","datatype","display.name") %in% ls())



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
for(i in 1:nrow(lod.peaks)) {
  
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



# Save the output of the find_peaks function
saveRDS(lod.peaks,file = paste0(prefix,'_lod_peaks_',threshold,'.rds'))


