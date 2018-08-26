#######################################################################################################################
# This script generates dataset.* for QTL viewer the attie metabolite data (cecum, liver, and plasma).
#
# Input:
#    1.) Path to .RData file generated from qtl2_gather_scan1_input_data.R
#    2.) Path to lod peaks .rds file generated from qtl2_findpeaks_scan1blup.R
#    3.) Datatype as character. Should be 'phenotype' as required by QTL viewer.
#    4.) Name for the dataset.* list.
#
# Output:
#    1.) A .RData file containing the dataset.*, genoprobs, markers, map, and kinship matrix
#
# Author: Duy Pham
# E-mail: duy.pham@jax.org
# Date: July 31, 2018
######################################################################################################################



### Read in the required data
#    1.) Path to .RData file generated from qtl2_gather_scan1_input_data.R
#    2.) Path to lod peaks .rds file generated from qtl2_findpeaks_scan1blup.R
#    3.) Datatype as character. Should be 'phenotype' as required by QTL viewer.
#    4.) Name for the dataset.* list.
load("~/Desktop/Attie/Metabolites/Cecum_Metabolites/attie_cecum_metabolite_qtl2_input.RData")
lod.peaks <- readRDS("~/Desktop/Attie/Metabolites/Cecum_Metabolites/attie_cecum_metabolite_rZ_lodpeaks_6.rds")
datatype <- 'phenotype'
dataset_name <- 'dataset.cecum.lipids'



### Formatiting the covar.factors 
covar.factors <- data.frame(column.names = covar.factors)
covar.factors$display.name <- gsub("([A-Z])([A-Z])", "\\1\\2 ", covar.factors$column.names)
covar.factors$display.name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", covar.factors$display.name, perl=TRUE)



### Changing raw dataframe to matrix
raw <- as.matrix(raw)



### Creating the formatted dataset as required by QTL Viewer
assign(dataset_name, list(annots = pheno.dict,
                          covar = covar,
                          covar.factors = covar.factors,
                          datatype = datatype,
                          display.name = display.name,
                          lod.peaks = list(additive = lod.peaks),
                          pheno = rankz,
                          norm = norm,
                          raw = raw,
                          samples = samples))



### Removing the data that is not stored in the dataset list
rm(dataset_name, covar, covar.factors, datatype, display.name, lod.peaks, norm, pheno.dict, rankz, raw, samples)
