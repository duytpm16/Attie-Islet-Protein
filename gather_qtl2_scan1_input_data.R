################################################################################
# Gather together the Attie liver metabolite, lipid, covariate and genoprobs
# data. Create objects taht can be fed into qtl2 for mapping including:
# pheno: data.frame containing covariates and phenotypes.
# pheno.descr: data.frame with pheotype descriptions.
# genoprobs: qtl2 style genoprobs object, one per chromosome.
# K: list containing kinship matrices, one per chromosome.
# map: list of marker positions, one per chromosome.
#
# We will create two Rdata files: one each for metabolites and lipids.
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
expr <- readRDS("attie_islet_protein_normalized.rds")      # Phenotype data
pheno.dict <- read.delim("~/Desktop/pQTL_Project/data/islet_proteins_pheno_dict.txt")  # Phenotype dictionary
genoprobs <- readRDS("~/Desktop/pQTL_Project/data/attie_DO500_genoprobs_20180303.rds") # Genoprobs data
markers <- readRDS("~/Desktop/pQTL_Project/data/marker_grid_0.02cM_plus.rds")          # Marker data
samples <- readRDS("attie_samples_annot.rds")              # Samples annotation data



### Checking if colnames in protein abundance and sample data matches the data dictionary (pheno.dict)
stopifnot(sum(ncol(samples),ncol(expr)) == nrow(pheno.dict))
stopifnot(c(colnames(samples),colnames(expr)) == pheno.dict$R_name)



### Subset genoprobes by mouse samples in normalized data.
mouse.samples  = sort(intersect(rownames(expr), rownames(genoprobs[[1]])))
for(i in 1:length(genoprobs)) {
  genoprobs[[i]] = genoprobs[[i]][mouse.samples,,]
}



### Check to see if the number of mouse samples in the protein abundance and genoprob data match
stopifnot(nrow(expr) == nrow(genoprobs[[1]]))



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
save(expr, pheno.dict, genoprobs, K, map, markers, covar, covar.factors, samples, file =  "attie_islet_proteins_qtl2_input.Rdata")


