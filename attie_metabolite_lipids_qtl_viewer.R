### Read in the required data
load("~/Desktop/Attie/Lipids/Cecum_Lipids/attie_cecum_lipid_qtl2_input.Rdata")
lod.peaks <- readRDS("~/Desktop/Attie/Lipids/Cecum_Lipids/attie_cecum_lipid_rZ_lod_peaks_6.rds")
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
