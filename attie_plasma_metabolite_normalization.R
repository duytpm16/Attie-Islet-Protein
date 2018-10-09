### This script takes in the raw plasma metabolite data file and normalizes 
#     the data by imputing missing values using pca and comBat normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
### Input:
#     1.) Prefix name to save the output files below.
#     2.) raw Attie plasma metabolite file: FilterdbyR_DOPlasmaMetbolites_BatchandCovariatesAppended.txt
#     3.) Attie sample annotation file: attie_DO_sample_annot.txt
#     4.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#
### Output:
#     1.) Filtered raw data file as rds file
#     2.) Normalized plasma metabolite levels by pca and comBat normaliziation method as a matrix in  rds file
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


### Options
options(stringsAsFactors = FALSE)



### Read in the raw plasma metabolite data, sample annotation, and Chr M and Y files.
#     plasma_metbaolite_raw: 380 x 342
#     samples: 500 x 7
#     chr: 498 x 5
prefix <- 'attie_plasma_metabolite'
raw <- read.delim("~/Desktop/Attie/Metabolites/Plasma_Metabolites/FilterdbyR_DOPlasmaMetbolites_BatchandCovariatesAppended.txt")
samples <- read.table("~/Desktop/Attie/attie_DO_sample_annot.txt", header = TRUE ,sep = "\t")
chr_m_y <- read.csv("~/Desktop/Attie/attie_sample_info_ChrM_Y.csv") 


### Variables to store the data
raw_file <- paste0(prefix,"_filtered_raw.rds")
norm_file <- paste0(prefix,"_normalized.rds")
norm_rz_file <- paste0(prefix,"_rZ_normalized.rds")
samples_file <-  paste0(prefix, "_samples_annot.rds")



### Fixing the name of two columns in the samples dataframe to match a data dictionary that will be used later in other scripts
colnames(samples)[grep('wave',colnames(samples), ignore.case = TRUE)] <- 'DOwave'
colnames(samples)[grep('batch',colnames(samples), ignore.case = TRUE)] <- 'batch'
colnames(raw)[grep('\\wave$', colnames(raw), ignore.case = TRUE)] <- 'DOwave'
colnames(raw)[grep('\\batch$', colnames(raw), ignore.case = TRUE)] <- 'batch'

raw$Mouse.ID <- gsub('-', '', raw$Mouse.ID)
samples$Mouse.ID <- gsub('-', '', samples$Mouse.ID)
chr_m_y$Mouse.ID <- gsub('-', '', chr_m_y$Mouse.ID)
colnames(samples) <- gsub('_','.',colnames(samples))



### Preparing samples dataframe
#   First, merge columns that are not metabolite to samples data.frame.
#   Next merge unique columns in chr_m_y dataframe to samples dataframe.
#       samples: 
samples <- merge(samples, raw[,c("Mouse.ID","DOwave","birthdate","sex","sac.date","coat.color","batch")], 
                 by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
samples <- samples[,-grep('\\.x$',colnames(samples))]
colnames(samples) <- gsub('\\.y$','',colnames(samples))
colnames(samples) <- gsub('_','.',colnames(samples))
colnames(samples)[grep('Mouse.ID',colnames(samples), ignore.case = TRUE)] <- 'mouse.id'



### Removing column that are not metabolite data
#     Inital dimensions: 380 x 342
#     After dimensions: 380 x 336 / 335 without Mouse.ID
raw <- raw[grep('DO', raw$Mouse.ID),]
rownames(raw) <- raw$Mouse.ID
raw <- raw[,!(colnames(raw) %in% c("Mouse.ID","DOwave","birthdate","sex","sac.date","coat.color","batch"))]



### Remove samples with more than 25% missing data.
#     Initial Dimensions: 380 x 336
#     After Dimensions: 368 x 336
prop.missing = rowMeans(is.na(raw))
sum(prop.missing > 0.25)
rownames(raw)[prop.missing > 0.25]

keep = which(prop.missing < 0.25)
raw = raw[keep,]
samples = samples[keep,]


### Log transformation of the protien abundance
data.log = log(raw[,!(colnames(raw) %in% 'Mouse.ID')])



### Set up batch and model for comBat
samples$sex  = factor(samples$sex)
samples$DOwave = factor(samples$DOwave)
mod = model.matrix(~sex, data = samples)
batch = samples$batch


### If there are no missing data, just combat normalize.
if(sum(is.na(raw) > 0)){
  
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
}else{
  data.log = ComBat(dat = t(data.log), batch = batch, mod = mod)
  data.log = t(data.log)
}



### Rank Z of normalized data.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

data.rz = data.log

for(i in 1:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}



### Saving the data to current working directory
saveRDS(raw, raw_file)
saveRDS(data.log, norm_file)
saveRDS(data.rz, norm_rz_file)
saveRDS(samples, samples_file)


