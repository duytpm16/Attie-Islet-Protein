####################################################################################################################################
#  This script takes in the raw cecum metabolite data file and normalizes 
#     the data by imputing missing values using pca and comBat normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
# 
#  Input:
#     1.) Prefix name to save the output files below.
#     2.) Raw Attie cecum metabolite file: FilterdbyR_DOCecumMetbolites_BatchandCovariatesAppended.txt
#     3.) Attie sample annotation file: attie_DO_sample_annot.txt
#     4.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#
#
#  Output:
#     1.) Filtered raw data file 
#     2.) Normalized cecum metabolite levels by pca and comBat normaliziation method as a matrix
#     3.) Normalized ranked z data 
#     4.) New Sample annotation dataframe 
#
#
#  Author: Duy Pham
#  Date:   July 10, 2018
#  E-mail: duy.pham@jax.org
####################################################################################################################################

### Load required libraries
options(stringsAsFactors = FALSE)
library(sva)
library(pcaMethods)









### Read in the raw cecum metabolite data, sample annotation, and Chr M and Y files.
#     Cecum metabolites raw: 366 x 507
#     samples: 500 x 7
#     chr: 498 x 5
raw <- read.delim("~/Desktop/Attie Mass Spectrometry/Metabolites/Cecum/FilterdbyR_DOCecumMetbolites_BatchandCovariatesAppended.txt")
samples <- read.table("~/Desktop/Attie Mass Spectrometry/Sample Info/attie_DO_sample_annot.txt", header = TRUE ,sep = "\t")
chr_m_y <- read.csv("~/Desktop/Attie Mass Spectrometry/Sample Info/attie_sample_info_ChrM_Y.csv") 











### Fixing the name of two columns in the samples dataframe to match a data dictionary that will be used later in other scripts
colnames(samples)[grep('wave',colnames(samples), ignore.case = TRUE)] <- 'DOwave'
colnames(samples)[grep('batch',colnames(samples), ignore.case = TRUE)] <- 'batch'
colnames(raw)[grep('wave', colnames(raw), ignore.case = TRUE)] <- 'DOwave'
colnames(raw)[grep('batch', colnames(raw), ignore.case = TRUE)] <- 'batch'

raw$Mouse.ID <- gsub('-', '', raw$Mouse.ID)
samples$Mouse.ID <- gsub('-', '', samples$Mouse.ID)
chr_m_y$Mouse.ID <- gsub('-', '', chr_m_y$Mouse.ID)
colnames(samples) <- gsub('_','.',colnames(samples))









### Preparing samples dataframe
#   First, merge columns that are not metabolite to samples data.frame.
#   Next merge unique columns in chr_m_y dataframe to samples dataframe.
#       samples: 366 x 11
samples <- merge(samples, raw[,c("Mouse.ID","DOwave","birthdate","sex","sac.date","coat.color","batch")], 
                 by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
samples <- samples[,-grep('\\.y$',colnames(samples))]
colnames(samples) <- gsub('\\.x$','',colnames(samples))
colnames(samples) <- gsub('_','.',colnames(samples))
colnames(samples)[grep('Mouse.ID',colnames(samples), ignore.case = TRUE)] <- 'mouse.id'







### Removing columns that are not phenotypes
#     Remove non-phenotype columns
#       Initial dimensions: 366 x 507
#       After dimensions: 366 x 500
raw <- raw[grep('DO', raw$Mouse.ID),]
rownames(raw) <- raw$Mouse.ID
raw <- raw[,!(colnames(raw) %in% c("Mouse.ID","DOwave","birthdate","sex","sac.date","coat.color","batch"))]







### 0 duplicates
sum(duplicated(rownames(raw)))








### Remove samples with more than 25% missing data.
#     Initial dimension: 366 x 500
#     After dimensiom: 364 x 500
prop.missing = rowMeans(is.na(raw))
sum(prop.missing > 0.25)
rownames(raw)[prop.missing > 0.25]


keep = which(prop.missing < 0.25)
raw = raw[keep,]
samples = samples[keep,]
rownames(samples) <- rownames(raw)









### Log transformation of the cecum metabolite
data.log = log(raw)









### Set up batch and model for comBat
samples$sex  = factor(samples$sex)
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









### Covariates
covar <- model.matrix(~ sex + DOwave + batch, data = samples)[,-1,drop = FALSE]

covar.info <- data.frame(sample.column = c('sex', 'DOwave', 'batch'),
                         covar.column  = c('sexM', 'DOwave', 'batch'),
                         display.name  = c('Sex', 'DO wave', 'Batch'),
                         interactive   = c(TRUE, FALSE, FALSE),
                         primary       = c(TRUE, FALSE, FALSE),
                         lod.peaks     = c('sex_int', NA, NA))









### QTL viewer format
dataset.cecum.metabolites <- list(annot.phenotype = data.frame(),
                                  annot.samples   = as_tibble(samples),
                                  covar.matrix    = covar,
                                  covar.info      = as_tibble(covar.info),
                                  data            = list(norm = data.log,
                                                    raw  = raw,
                                                    rz   = data.rz),
                                  datatype        = 'phenotype',
                                  display.name    = 'Attie Cecum Metabolites',
                                  lod.peaks       = list())









### Save
rm(list = ls()[!grepl('dataset[.]', ls())])
save(dataset.cecum.metabolites, file = '~/Desktop/Attie Mass Spectrometry/Metabolites/Cecum/attie_cecum_metabolite_viewer.Rdata')



