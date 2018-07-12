### This script gathers the chunked QTL files outputted from the scan1 function 
#      and saves chunk's QTL LOD as one large matrix
#
### Input:
#      1.) input.dir:     directory to where the QTL chunks are stored
#      2.) output.prefix: prefix of the file name to save the large qtl matrix
#      3.) pheno.dict:    path to data dictionary
#
### Output:
#      1.) qtl.mat: saves large QTL LOD of all chunks to current directory.
#
### Author: Duy Pham and Dan Gatti
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org, dan.gatti@jax.org
#######################################################################################################################



### Install required packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl","dplyr"))
# library(devtools)
# install_github("rqtl/qtl2")
# install_github("rqtl/qtlconvert")

### Load required libraries
library(qtl2)
library(qtl2convert)

### Options
options(stringsAsFactors = F)



### Variables to change
#     input.dir:     Full path to directory to chunks of qtl, no '~'
#     output.prefix: Prefix to name of file to save all qtl data
#     pheno.dict:    Data dictionary file. Just to double check if we have all the proteins.
input.dir <- getwd() 
output.prefix <- "attie_islet_protein"
pheno.dict <- read.delim("islet_proteins_pheno_dict.txt")



### Get the QTL files, which may contain more than one result.
qtl.files <- dir(path = input.dir, pattern = "_QTL.rds$", full.names = T)



### Verify that we have all of the QTL files.
chunk.num = gsub(input.dir, "", qtl.files)
chunk.num = gsub(".*chunk_|_QTL.rds$", "", chunk.num)
chunk.num = sort(as.numeric(chunk.num))



# Check to see if we have all the chunks
sd <- setdiff(1:max(chunk.num), chunk.num)
if(length(sd) > 0) {
  
  stop(paste("Some output files are missing:", paste(sd, collapse = ", ")))
  
} 



# Goes through each chunk file and save all QTLs to one matrix, qtl.mat

qtl.mat = NULL

for(i in 1:length(qtl.files)) {
  
  print(paste("Chunk", i, "of", length(qtl.files)))
  
  # Read in one of the chunk file
  qtl = readRDS(qtl.files[i])
  
  # Add these QTLs to the large QTL matrix.
  if(is.null(qtl.mat)) {
      qtl.mat = qtl
  } else {
      qtl.mat = cbind(qtl.mat, qtl)
  } 
  
}



### Convert large QTL LOD matrix to class scan1 and matrix
class(qtl.mat) = c("scan1", "matrix")



### Check to see if we have all the proteins
protein.names <- as.character(pheno.dict[pheno.dict$is_pheno, 'data_name'])
stopifnot(protein.names == colnames(qtl.mat))



### Save the large QTL LOD matrix to current directory
saveRDS(qtl.mat, file = paste0(output.prefix, "_all_lod_qtl.rds"))



### Remove unnecessary data
rm(chunk.num, fig.dir, i, input.dir, output.prefix, qtl.files, sd, qtl, protein.names)
