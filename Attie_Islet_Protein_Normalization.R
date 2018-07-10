### This script takes in the raw protein abundance data file and normalizes 
#     the data by imputing missing values using pca and comBat normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
### Input:
#     1.) raw Attie islet protein abundance file: DO_islet_proteomics_non_normalized.txt
#     2.) Attie sample annotation file: attie_DO_sample_annot.txt
#     3.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#
### Output:
#     Two rds file containing:
#         1.) Normalized protein abundance levels by pca and comBat normaliziation method as a matrix: attie_islet_protein_normalized
#         2.) Sample annotation dataframe: samples.
#
### Author: Duy Pham, most of the codes were taken from Dan Gatti's script
### Date:   July 10, 2018
#####################################################################



### Install required packages
# source("https://bioconductor.org/biocLite.R")
# bioclite(c('sva','pcaMethods'))
# install.packages('tidyverse')

### Load required libraries
library(sva)
library(pcaMethods)
library(tidyverse)

### Options
options(stringsAsFactors = FALSE)



### Read in the raw islet proteomic data, sample annotation, and Chr M and Y files.
#     islet_protein_raw: 439 x 8239
#     samples: 500 x 7
#     chr: 498 x 5
islet_protein_raw <- read.table("DO_islet_proteomics_non_normalized.txt", sep = '\t', header = TRUE)
samples <- read.table("attie_DO_sample_annot.txt", header = TRUE ,sep = "\t")
chr_m_y <- read.csv("attie_sample_info_ChrM_Y.csv") 



### Preparing samples dataframe
#   First, merge columns that do not contain protein abundance to samples data.frame.
#   Next merge unique columns in chr_m_y dataframe to samples dataframe.
#       samples: 375 x 10 (Reduced to 375 because control ('Std) mice were removed)
samples <- merge(samples, islet_protein_raw[,c("Mouse.ID","Injection_order","Plate_number","Batch")], by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
samples$Mouse.ID <- gsub('-','',samples$Mouse.ID)



### Keep protein abundance columns that have less than 50% NAs across samples.
#       islet_protein_raw: 439 x 5434
islet_protein_raw <- islet_protein_raw[,!(colnames(islet_protein_raw) %in% 
                                            c("Injection_order","Plate_number","Batch"))]
islet_protein_raw <-  islet_protein_raw[,colSums(is.na(islet_protein_raw)) < .50 * nrow(islet_protein_raw)]
islet_protein_raw$Mouse.ID <- gsub ('-', '', islet_protein_raw$Mouse.ID)



### Remove control samples.
#       Number of controls: 64 with some duplicates
#       Number of DOs: 375 with some duplicates (DO-174, DO-374)
ctrl <- islet_protein_raw[grep("Std", islet_protein_raw$Mouse.ID),]
islet_protein_raw <- islet_protein_raw[grep('DO', islet_protein_raw$Mouse.ID),]



### Log transformation of the protien abundance
data.log = log(islet_protein_raw[,!(colnames(islet_protein_raw) %in% 'Mouse.ID')])



### Set up batch and model for comBat
samples$sex  = factor(samples$sex)
samples$wave = factor(samples$wave)
mod = model.matrix(~sex, data = samples)
batch = samples$Batch



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
    
    data.log = data.cb
    break
    
  }
})



### Find the duplicated samples and keep the one with fewer NAs in original data
dupl = which(duplicated(islet_protein_raw$Mouse.ID))
prop.missing = rowMeans(is.na(islet_protein_raw))
unique.samples = unique(islet_protein_raw$Mouse.ID)
keep = rep(FALSE, nrow(islet_protein_raw))

for(i in 1:length(unique.samples)){
  
  sample = unique.samples[i]
  wh = which(islet_protein_raw$Mouse.ID == sample)
  wh = wh[which.min(prop.missing[wh])]
  keep[wh] = TRUE

} 

data.log <- data.log[keep,]
data.log <- sapply(data.log, as.numeric)
rownames(data.log) <- islet_protein_raw[keep,'Mouse.ID']



### Removing duplicates in the samples dataframe and removing '-' in Mouse ID
samples = samples[match(rownames(data.log), samples$Mouse.ID),]
rownames(samples) <- samples$Mouse.ID



### Saving the data to current working directory
saveRDS(data.log, "attie_islet_protein_normalized.rds")
saveRDS(samples, "attie_samples_annot.rds")



### Remove other data
rm(chr_m_y,ctrl,data.cb,mod,pc.data,batch,chg,dupl,i,iter,keep,miss,prop.missing,sample,unique.samples, data.compl, islet_protein_raw, wh)


