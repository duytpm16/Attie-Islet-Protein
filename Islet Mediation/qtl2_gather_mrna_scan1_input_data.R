### This script gathers the Attie expression data (raw filtered, normalized, and rankz normalized), 
#    marker maps, and genoprobs files. 
#    Addtionally, this script create objects will be inputted into qtl2 for mapping and QTL Viewer including:
#         covar: data matrix containing covariates
#         K: list containing kinship matrices, one per chromosome.
#         map: list of marker positions, one per chromosome.
#
### Input:
#     1.) filtered raw  data
#     2.) normalized islet mrna data by upper quantile (.75)
#     3.) rankz of the normalized islet data
#     4.) genoprobs data
#     5.) markers dataframe
#     6.) samples annotation data
#     7.) datatype: Should be "protein" or "mRNA"
#     8.) display.name: Name to display on QTL viewer
#     9.) Rdata file name to store all the data 
#
### Output:
#     1.) Rdata file of all data above.
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
prefix <- "~/Desktop/Islet Mediation/mRNA/attie_islet_mrna_284"                        # Prefix of the new raw, normalized, rank Z normalized, and samples rds file
genoprobs <- readRDS("~/Desktop/Attie/attie_DO500_genoprobs_20180303.rds")             # Genoprobs data
markers <- readRDS("~/Desktop/Attie/marker_grid_0.02cM_plus.rds")                      # Marker data
datatype <- "mRNA"                                                                     # Datatype for QTL viewer
display.name <- "Attie Islet mRNA 284"                                                 # Display name for QTL viewer



### Load in the other data.
raw <- readRDS(paste0(prefix,"_filtered_raw.rds"))         # Raw protein levels
norm <- readRDS(paste0(prefix,"_normalized.rds"))          # Normalized protein levels by comBat
rankz <- readRDS(paste0(prefix, "_rZ_normalized.rds"))     # RankZ normalized protein levels
samples <- readRDS(paste0(prefix, "_samples_annot.rds"))   # Samples annotation data
gathered_data <- paste0(prefix,"_qtl2_input.RData")        # Rdata file name to store all data




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




### Set up covariates. 
#     There may be NA's in the covar.names, which will be removed.
covar.factors <- c('sex','DOwave')




### Check to see if the covar.factors matches the ones in samples
stopifnot(covar.factors %in% colnames(samples))



### Covar matrix
f = as.formula(paste("~", paste(covar.factors, collapse = "+")))
covar <- model.matrix(f, data = samples)[,-1]



### Calculate kinship matrix.
K = calc_kinship(probs = genoprobs, type = "loco", quiet = FALSE, cores = 4)


### Removing some data
rm(f, i, mouse.samples, prefix)

### Save to *.Rdata file.
save(norm, genoprobs, K, map, markers, covar, covar.factors, samples, raw, rankz, datatype, display.name,
     file =  gathered_data)


