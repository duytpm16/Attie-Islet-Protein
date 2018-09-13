### This script gathers the Attie expression data (raw filtered, normalized, and rankz normalized),
#    samples dataframe, marker maps, genoprobs, and phenotype dictionary files. 
#    Addtionally, this script will create 3 new objects:
#         covar: data matrix containing covariates
#         K: list containing kinship matrices, one per chromosome.
#         map: list of marker positions, one per chromosome.
#    All of these objects will be saved as a .RData file used for the scan1 function in the qtl2 package.
#
### Input:
#     1.) Prefix name that was given in the normalization script to get the 
#             filtered raw, normalized, rankz, and samples dataframe.
#         This prefix name will also be used to save the output .RData file as prefix + '_qtl2_input.RData'
#     2.) phenotype dictionary file
#     3.) genoprobs data file
#     4.) markers dataframe file
#     5.) datatype: Should be "protein" or "mRNA"
#     6.) display.name: Name to display on QTL viewer
#
### Output:
#     1.) Rdata file of all data above
# 
### Author: Duy Pham and Dan Gatti
### Date:   July 11, 2018
### E-mail: duy.pham@jax.org,dan.gatti@jax.org
#
################################################################################



### Install required packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl","dplyr"))
# library(devtools)
# install_github("rqtl/qtl2")
# install_github("rqtl/qtlconvert")

### Load required libraries
library(qtl2)
library(qtl2convert)
library(tidyverse)

### Options
options(stringsAsFactors = F)


  
### Variables to change 
prefix <- "attie_plasma_metabolite"                            # Prefix of the new raw, normalized, rank Z normalized, and samples rds file
pheno.dict <- read.delim("plasma_metabolites_pheno_dict.txt")  # Phenotype data dictionary
genoprobs <- readRDS("attie_DO500_genoprobs_20180303.rds")     # Genoprobs data
markers <- readRDS("marker_grid_0.02cM_plus.rds")              # Marker data
datatype <- "mRNA"                                             # Datatype for QTL viewer
display.name <- "Attie Plasma Metabolite"                      # Display name for QTL viewer



### Creating name to save the output .RData file
gathered_data <- paste0(prefix,"_qtl2_input.RData")        # Rdata file name to store all data



### Load in the other data.
raw <- readRDS(paste0(prefix,"_filtered_raw.rds"))         # Raw protein levels
norm <- readRDS(paste0(prefix,"_normalized.rds"))          # Normalized protein levels by comBat
rankz <- readRDS(paste0(prefix, "_rZ_normalized.rds"))     # RankZ normalized protein levels
samples <- readRDS(paste0(prefix, "_samples_annot.rds"))   # Samples annotation data




# Preparing pheno.dict names
for(i in 1:ncol(pheno.dict)){
  if(is.logical(pheno.dict[,i])){
    break
  }
  pheno.dict[!pheno.dict$is_pheno,i] <- gsub('_','.', pheno.dict[!pheno.dict$is_pheno,i])
}
pheno.dict[pheno.dict == 'Mouse.ID'] <- 'mouse.id'
samples <- samples[,match(colnames(samples), pheno.dict[!pheno.dict$is_pheno, 'data_name'])]



### Checking if colnames in protein abundance and sample data matches the data dictionary (pheno.dict)
stopifnot(sum(ncol(samples),ncol(norm)) == nrow(pheno.dict))
stopifnot(c(colnames(samples),colnames(norm)) == pheno.dict$data_name)



### Subset genoprobes by mouse samples in normalized data.
mouse.samples  = sort(intersect(rownames(norm), rownames(genoprobs[[1]])))
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][mouse.samples,,]
}



### Check to see if the number of mouse samples in the protein abundance and genoprob data match
stopifnot(nrow(norm) == nrow(genoprobs[[1]]))



### Make marker map.
map = map_df_to_list(map = markers[,c('marker','chr','pos','cM')], pos_column = "pos")



### Check to see if number of markers per chromosome in genoprobes and map matches
stopifnot(sum(sapply(genoprobs, dim)[3,]) == sum(sapply(map, length)))



### Create factors for covariates that we want to map as factors.
samples$sex    = factor(samples$sex)
samples$DOwave = factor(samples$DOwave)
samples$batch  = factor(samples$batch)



### Set up covariates. 
#     There may be NA's in the covar.names, which will be removed.
covar.factors <- strsplit(pheno.dict$covar_list, ":")
covar.factors <- unique(unlist(covar.factors))
covar.factors <- covar.factors[!is.na(covar.factors)]



### Check to see if the covar.factors matches the ones in samples
stopifnot(covar.factors %in% colnames(samples))



### Covar matrix
f = as.formula(paste("~", paste(covar.factors, collapse = "+")))
covar <- model.matrix(f, data = samples)[,-1]



### Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)



### Removing some data
rm(f, i, mouse.samples, prefix)



### Save to *.Rdata file
save(norm, pheno.dict, K, map, markers, covar, covar.factors, samples, raw, rankz, datatype, display.name, genoprobs,
     file =  gathered_data)


