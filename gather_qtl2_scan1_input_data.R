### This script gathers the Attie islet protein level data (raw filtered, normalized, and rankz normalized), 
#    marker maps, genoprobs, and phenotype dictionary files. 
#    Addtionally, this script create objects will be inputted into qtl2 for mapping and QTL Viewer including:
#         covar: data matrix containing covariates
#         K: list containing kinship matrices, one per chromosome.
#         map: list of marker positions, one per chromosome.
#
### Input:
#     1.) filtered raw islet protein data
#     2.) normalized islet protein data by pca and comBat
#     3.) rankz of the normalized islet protein data
#     4.) phenotype dictionary
#     5.) genoprobs data
#     6.) markers dataframe
#     7.) samples annotation data
#     8.) datatype: Should be "protein"
#     9.) display.name: Name to display on QTL viewer
#
### Output:
#     1.) Rdata file of all data above.
# 
# Daniel Gatti
# dan.gatti@jax.org
# July 17, 2017
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



### Load in the data.
raw <- readRDS("attie_islet_protein_filtered_raw.rds")     # Raw protein levels
norm <- readRDS("attie_islet_protein_normalized.rds")      # Normalized protein levels by comBat
rankz <- readRDS("attie_islet_protein_rZ_normalized.rds")  # RankZ normalized protein levels
pheno.dict <- read.delim("islet_proteins_pheno_dict.txt")  # Islet phenotype data dictionary
genoprobs <- readRDS("attie_DO500_genoprobs_20180303.rds") # Genoprobs data
markers <- readRDS("marker_grid_0.02cM_plus.rds")          # Marker data
samples <- readRDS("attie_samples_annot.rds")              # Samples annotation data
datatype <- "protein"                                      # Datatype for QTL viewer
display.name <- "Attie Islet Protein"                      # Display name for QTL viewer


### Checking if colnames in protein abundance and sample data matches the data dictionary (pheno.dict)
stopifnot(sum(ncol(samples),ncol(norm)) == nrow(pheno.dict))
stopifnot(c(colnames(samples),colnames(norm)) == pheno.dict$R_name)



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



### Remove unnecessary data
rm(f, i, mouse.samples)



### Save to *.Rdata file.
save(norm, pheno.dict, genoprobs, K, map, markers, covar, covar.factors, samples, raw, rankz, datatype, display.name,
     file =  "attie_islet_proteins_qtl2_input.Rdata")


