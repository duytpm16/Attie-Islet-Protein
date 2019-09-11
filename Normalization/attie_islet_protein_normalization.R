##########################################################################################################################
#
#  This script takes in the raw protein abundance data file and normalizes 
#     the data by imputing missing values using pca and comBat for normalization.
#     Additionally, a samples dataframe is created by reading two additional files below.
#
#
#  Input:
#     1.) Raw Attie islet protein abundance file: DO_islet_proteomics_non_normalized.txt
#     2.) Attie sample annotation file: attie_DO_sample_annot.txt
#     3.) Sample's Chr M and Y info file: "attie_sample_info_ChrM_Y.csv"
#     4.) Protein annotation file
#     5.) Genotype probabilities in qtl2 format
#     6.) 69,0005 marker grid
#
#  Output:
#     1.) .Rdata file in qtl2 viewer format
#
#
#
#  Author: Duy Pham
#  Date:   July 10, 2018
#  E-mail: duy.pham@jax.org
##########################################################################################################################


### Options
options(stringsAsFactors = FALSE)

### Load required libraries
library(sva)
library(qtl2)
library(dplyr)
library(pcaMethods)
library(qtl2convert)








### Variables to change
#     1.) Raw protein data
#     2.) Sample annotations
#     3.) Chr M and Y files.
#     4.) Annotation for the protein Uniprot ID
#     5.) Genotype probabilities in qtl2 format
#     6.) 69,0005 marker grid
#         
#     Variables          Raw dimension
#     ------------------ -----------------
#     islet_protein_raw: 439 x 8239
#     samples:           500 x 7
#     chr:               498 x 5
#     annots:            8,235 x 11
#     genoprobs:         500 x 8 x 69,005
#     markers:           500 x 7
raw       <- read.table("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Original Files/DO_islet_proteomics_non_normalized.txt", sep = '\t', header = TRUE)
samples   <- read.table("~/Desktop/Attie Mass Spectrometry/Sampe Info/attie_DO_sample_annot.txt", header = TRUE ,sep = "\t")
chr_m_y   <- read.csv("~/Desktop/Attie Mass Spectrometry/Sampe Info/attie_sample_info_ChrM_Y.csv") 
annots    <- readRDS("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Original Files/annotated_uniprotID_from_Coon_lab.rds")
genoprobs <- readRDS("~/Desktop/Attie Mass Spectrometry/Genotype/attie_DO500_genoprobs_20180303.rds")
markers   <- readRDS("~/Desktop/Attie Mass Spectrometry/Genotype/marker_grid_0.02cM_plus.rds")








### Fixing the name of two columns in samples dataframe
#   Remove hyphen in DO mouse ID
colnames(samples)[grep('wave', colnames(samples))] <- 'DOwave'
colnames(raw)[grep('Batch', colnames(raw))]        <- 'batch'

samples$Mouse.ID  <- gsub('-', '', samples$Mouse.ID)
raw$Mouse.ID      <- gsub('-', '', raw$Mouse.ID)
chr_m_y$Mouse.ID  <- gsub('-', '', chr_m_y$Mouse.ID)









### Preparing samples dataframe
#   First, merge columns that do not contain protein abundance to samples data.frame.
#       This includes the mouse.id, injection order, plate number, and batch.
#   Next merge columns in chr_m_y dataframe to samples dataframe.
#       samples: 375 x 10 (Reduced to 375 because control ('Std) mice were removed)
samples <- merge(samples, raw[,c("Mouse.ID","Injection_order","Plate_number","batch")], by = "Mouse.ID")
samples <- merge(samples, chr_m_y[,c('Mouse.ID','generation','chrM','chrY')], by = "Mouse.ID")
colnames(samples) <- gsub('_','.',colnames(samples))










### Remove control samples.
#       Number of controls: 64 with some duplicates
#       Number of DOs: 375 with two duplicates (DO-174, DO-374)
ctrl <- raw[grep("Std", raw$Mouse.ID),]
raw  <- raw[grep('DO', raw$Mouse.ID),]










### Find the duplicated samples and keep the one with fewer NAs in original data (Taken from D. Gatti)
#     raw: 373 x 8,235
#     sample: 373 x 13
dupl <- which(duplicated(raw$Mouse.ID))
prop.missing   <- rowMeans(is.na(raw))
unique.samples <- unique(raw$Mouse.ID)
keep = rep(FALSE, nrow(raw))

for(i in 1:length(unique.samples)){
  
  sample <- unique.samples[i]
  wh <- which(raw$Mouse.ID == sample)
  wh <- wh[which.min(prop.missing[wh])]
  keep[wh] = TRUE
} 

raw     <- raw[keep,]
samples <- samples[keep,]




### Make row names the mouse ID
rownames(raw)     <- raw$Mouse.ID
rownames(samples) <- samples$Mouse.ID
stopifnot(rownames(raw) == rownames(samples))
colnames(samples)[grep('mouse.id', colnames(samples), ignore.case = TRUE)] <- 'mouse.id'










### Removing columns that are not proteins
#      raw: 373 x 8,235
#   Keep protein columns that have less than 50% NAs across samples.
#      raw: 373 x 5,415
raw <- raw[,!(colnames(raw) %in% c("Mouse.ID","Injection_order","Plate_number","batch"))]
raw <- raw[,colSums(is.na(raw)) < .50 * nrow(raw)]






### Creating annotations for the proteins
#    Removed proteins that were not in raw dataframe
#    Removed proteins with missing gene id and protein id
#    Removed proteins with multiple gene ids
#    Keep proteins part of the autosome, sex, and mitochondrial chromosomes
#
#    annots: 4,818 x 11
annots$Majority.protein.IDs <- gsub(';','_',annots$Majority.protein.IDs, fixed = TRUE)
annots$Majority.protein.IDs <- gsub('-','.',annots$Majority.protein.IDs, fixed = TRUE)
annots <- annots[annots$Majority.protein.IDs %in% colnames(raw),]
annots <- annots[!is.na(annots$gene_id) & !is.na(annots$protein_id), ]
annots <- annots[!grepl(';', annots$gene_id), ]

annots <- annots %>%
            filter(chr %in% c(1:19,'X','Y','MT')) %>%
            dplyr::rename(uniprot_id = Majority.protein.IDs) %>%
            arrange(protein_id) %>%
            `rownames<-`(.$protein_id)






### Keep proteins with annotations in raw dataframe
#    raw: 373 x 4,818
raw <- raw[,colnames(raw) %in% annots$uniprot_id]
colnames(raw) <- annots$protein_id[match(colnames(raw), annots$uniprot_id)]



### Log transformation of the protein abundance
#
#   data.log: 373 x 4,818
data.log = log(raw)











### Set up batch and model for comBat
samples$sex    <- factor(samples$sex)
mod   <- model.matrix(~sex, data = samples)
batch <- samples$batch



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

data.rz = apply(data.log, 2, rankZ)








### Setting of covariates, kinship, and list for qtl2
covar <- model.matrix(~ sex + DOwave + batch, samples)[,-1]
K     <- calc_kinship(genoprobs, type = 'loco', cores = 10)
map   <- map_df_to_list(markers, pos_column = 'pos')











### QTL Viewer dataset
dataset.islet.proteins <- list(annot.protein = annots,
                               covar.matrix  = covar,
                               covar.factors = data.frame(),
                               data          = list(norm = data.log,
                                                    raw  = raw,
                                                    rz   = data.rz),
                               datatype      = 'protein',
                               display.name  = 'DO Islet Proteins',
                               lod.peaks     = list(),
                               samples       = samples)
ensembl.version = 91




rm(list = ls()[!ls() %in% c('dataset.islet.proteins','K','genoprobs','map','markers','enseml.version')])


save.image(file = 'attie_islet_proteins_qtl_viewer_v1.RData')

