### This script takes in the raw liver lipids data file and normalizes 
#     the data by imputing missing values using pca and comBat normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
### Input:
#     1.) raw Attie liver lipids file: FilterdbyR_DOCecumMetbolites_BatchandCovariatesAppended.txt
#     2.) Attie sample annotation file: attie_DO_sample_annot.txt
#     3.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#     4.) Name to store the filtered raw data as rds file
#     5.) Name to store the normalized data as rds file
#     6.) Name to store the normalized rank z data as rds file
#     7.) Name to store the new sample annotation data as rds file
#
### Output:
#     1.) Filtered raw data file as rds file
#     2.) Normalized liver lipid levels by pca and comBat normaliziation method as a matrix in  rds file
#     3.) Normalized ranked z data as rds file
#     4.) New Sample annotation dataframe as rds file
#
### Author: Duy Pham, most of the codes were taken from Dan Gatti's script
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org
#####################################################################


### Variables to change
#      liver_lipid_raw: 419 x 1564
#      samples: 500 x 7
#      chr_m_y: 498 x 5
prefix <- '~/Desktop/Attie/Lipids/Liver_Lipids/attie_liver_lipid'
raw <- read.delim("~/Desktop/Attie/Lipids/Liver_Lipids/03_January_2018_DO_Liver_Lipidomics_Raw.txt")
samples <- read.delim('~/Desktop/Attie/attie_DO_sample_annot.txt')
chr_m_y <- read.csv("~/Desktop/Attie/attie_sample_info_ChrM_Y.csv")



### Variable names to store data
raw_file <- paste0(prefix,"_filtered_raw.rds")
norm_file <- paste0(prefix,"_normalized.rds")
norm_rz_file <- paste0(prefix,"_rZ_normalized.rds")
samples_file <-  paste0(prefix, "_samples_annot.rds")



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
#     First, merge columns that are not lipids to samples data.frame.
#     Next merge unique columns in chr_m_y dataframe to samples dataframe.
#         Dimension: 384 x 11
samples <- merge(samples, raw[,c("Mouse.ID","batch")], by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
colnames(samples) <- gsub('_','.',colnames(samples))
colnames(samples)[grep('Mouse.ID',colnames(samples), ignore.case = TRUE)] <- 'mouse.id'


### Removing controls from raw dataframe and removing non-phenotype columns (except Mouse.ID):
#     Initial dimensions: 439 x 1564
#     After removing controls: 384 x 1564
#     After removing non-phenotype columns: 384 x 1562
raw <- raw[grep('DO', raw$Mouse.ID),]
rownames(raw) <- raw$Mouse.ID
raw <- raw[,!(colnames(raw) %in% c("Mouse.ID","DOwave","batch"))]



### Swapping DO 343/373
DO343 <- raw["DO373",]
DO373 <- raw["DO343",]
raw["DO343",] <- DO343[1,]
raw["DO373",] <- DO373[1,]


### 0 NAs and duplicates 
sum(is.na(raw))
sum(duplicated(rownames(raw)))
rownames(samples) <- rownames(raw)


### Log transformation of the liver lipids
data.log = log(raw)



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



### Remove other data
rm(chr_m_y,data.cb,mod,pc.data,batch,chg,dupl,i,rankZ, DO343, DO373,
   iter,keep,miss,prop.missing,sample,unique.samples, data.compl, wh, norm_file, norm_rz_file, raw_file, samples_file, prefix)
