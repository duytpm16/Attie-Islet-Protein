### This script takes in the raw protein abundance data file and normalizes 
#     the data by imputing missing values using pca and comBat for normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
### Input:
#     1.) raw Attie islet protein abundance file: DO_islet_proteomics_non_normalized.txt
#     2.) Attie sample annotation file: attie_DO_sample_annot.txt
#     3.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#     4.) Name to store the filtered raw data as rds file
#     5.) Name to store the normalized data as rds file
#     6.) Name to store the normalized rank z data as rds file
#     7.) Name to store the new sample annotation data as rds file
#
### Output:
#     1.) Filtered raw data file as rds file
#     2.) Normalized protein abundance levels by pca and comBat normaliziation method as a matrix in  rds file
#     3.) Normalized ranked z data as rds file
#     4.) New Sample annotation dataframe as rds file
#
### Author: Duy Pham, most of the codes were taken from Dan Gatti's script
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org
#####################################################################



### Install required packages
# source("https://bioconductor.org/biocLite.R")
# bioclite(c('sva','pcaMethods'))
# install.packages('tidyverse')

### Load required libraries
library(sva)
library(pcaMethods)
library(dplyr)


### Options
options(stringsAsFactors = FALSE)



### Variables to change
#     1.) Prefix to describe that describes the data and will be used to get files in the pipeline
#     2.) sample annotation file for the data
#     3.) Chr M and Y files.
#     4.) RNAseq data    
#     islet_protein_raw: 439 x 8239
#     samples:           500 x 7
#     chr:               498 x 5
prefix  <- "attie_islet_protein_284"
raw     <- read.table("~/Desktop/Attie Final/Islet/Proteins/DO_islet_proteomics_non_normalized.txt", sep = '\t', header = TRUE)
samples <- read.table("~/Desktop/Attie Final/attie_DO_sample_annot.txt", header = TRUE ,sep = "\t")
chr_m_y <- read.csv("~/Desktop/Attie Final/attie_sample_info_ChrM_Y.csv") 
load("~/Desktop/Attie_DO378_eQTL_viewer_v6.Rdata")       # Cadillac: /projects/churchill-lab/data/Attie/DO_islet_RNA



### Variable names to store the data
raw_file     <- paste0(prefix,"_filtered_raw.rds")
norm_file    <- paste0(prefix,"_normalized.rds")
norm_rz_file <- paste0(prefix,"_rZ_normalized.rds")
samples_file <- paste0(prefix, "_samples_annot.rds")











### Fixing the name of two columns: 1 in the samples and 1 in raw dataframe to match a data dictionary that will be used later in other scripts
samples <- samples %>% dplyr::rename(DOwave = wave)
raw     <- raw %>% dplyr::rename(batch = Batch)
samples$Mouse.ID  <- gsub('-', '', samples$Mouse.ID)
raw$Mouse.ID      <- gsub('-', '', raw$Mouse.ID)
chr_m_y$Mouse.ID  <- gsub('-', '', chr_m_y$Mouse.ID)


### Preparing samples dataframe
#   First, merge columns that do not contain protein abundance to samples data.frame.
#   Next merge unique columns in chr_m_y dataframe to samples dataframe.
#       samples: 375 x 10 (Reduced to 375 because control ('Std) mice were removed)
samples <- merge(samples, raw[,c("Mouse.ID","Injection_order","Plate_number","batch")], by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
colnames(samples) <- gsub('_','.',colnames(samples))











### Remove control samples.
#       Number of controls: 64 with some duplicates
#       Number of DOs: 375 with some duplicates (DO-174, DO-374)
ctrl <- raw[grep("Std", raw$Mouse.ID),]
raw  <- raw[grep('DO',  raw$Mouse.ID),]













### Find the duplicated samples and keep the one with fewer NAs in original data
dupl = which(duplicated(raw$Mouse.ID))
prop.missing = rowMeans(is.na(raw))
unique.samples = unique(raw$Mouse.ID)
keep = rep(FALSE, nrow(raw))

for(i in 1:length(unique.samples)){
  
  sample = unique.samples[i]
  wh = which(raw$Mouse.ID == sample)
  wh = wh[which.min(prop.missing[wh])]
  keep[wh] = TRUE
} 

raw <- raw[keep,]
samples = samples[keep,]


### Make row names the mouse ID
rownames(raw)     <- raw$Mouse.ID
rownames(samples) <- samples$Mouse.ID
colnames(samples)[grep('mouse.id', colnames(samples), ignore.case = TRUE)] <- 'mouse.id'











### Find overlapping samples between rnaseq and protein data (284 samples)
overlap.samples <- intersect(rownames(raw), rownames(dataset.islet.rnaseq$raw))
raw <- raw[overlap.samples,]
samples <- samples[overlap.samples,]











### Removing columns that are not proteins
#   Keep protein columns that have less than 50% NAs across samples.
#   284 x 8239 -> 284 x 8235 -> 284 x 5469
#    
raw <- raw[,!(colnames(raw) %in% c("Mouse.ID","Injection_order","Plate_number","batch"))]
raw <- raw[,colSums(is.na(raw)) < .50 * nrow(raw)]












### Log transformation of the protein abundance
#
#    284 x 5469
data.log = log(raw)











### Set up batch and model for comBat
samples$sex  = factor(samples$sex)
samples$DOwave = factor(samples$DOwave)
mod = model.matrix(~sex, data = samples)
batch = samples$batch


### Impute missing data and combat normalization
chg = 1e6
iter = 1
repeat({
  
    print(paste("Iteration", iter))
    
    # Impute missing data using pca
    miss = which(is.na(data.log))
    print(paste(length(miss), "missing points."))
    
    pc.data = pca(data.log, method = "bpca", nPcs = 5)
    data.compl = completeObs(pc.data)
    
    # Batch adjust.
    # ComBat wants the data with variable in rows and samples in columns.
    data.cb = ComBat(dat = t(data.compl), batch = batch, mod = mod)
    data.cb = t(data.cb)
    
    # Calculate the change.
    chg = sum((data.compl[miss] - data.cb[miss])^2)
    print(paste("SS Change:", chg))
    
    # Put the missing data back in and impute again.
    if(chg > 5 & iter < 20) {
      
       data.cb[miss] = NA
       data.log = data.cb
       iter = iter + 1    
      
    }else{
      
       data.cb[miss] = NA
       data.log = data.cb
       break
      
    }
})




  
  
  
  
  
  
### Rank Z of normalized data.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


data.rz = data.log
data.rz = apply(data.rz, 2, rankZ)










### Saving the data to current working directory
saveRDS(raw, raw_file)
saveRDS(data.log, norm_file)
saveRDS(data.rz, norm_rz_file)
saveRDS(samples, samples_file)
