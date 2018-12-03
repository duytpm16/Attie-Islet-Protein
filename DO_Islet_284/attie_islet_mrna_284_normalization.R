#################################################################################################################
#  This script is used to re-normalize the attie islet mrna data based on samples that overlap with 
#     attie islet protein.
#
#  Input:
#    1.) .RData file that contains dataset.islet.proteins and dataset.islet.rnaseq (Attie_QTL_Viewer_V2.RData)
#  
#  Output:
#    1.) .rds raw data
#    2.) .rds normalized data by upper quantile
#    3.) .rds rankz data 
#    4.) .rds sample
#
#  Author: Duy Pham
#  E-mail: duy.pham@jax.org
#################################################################################################################

### Load in Attie QTL viewer data
load("~/Desktop/Attie/Attie_rZ_QTL_Viewer_V2.RData")



### Find overlapping samples between islet protein and mrna datasets
intersect.samples <- intersect(rownames(dataset.islet.rnaseq$expr), rownames(dataset.islet.proteins$expr))



### Extract raw and samples data based on overlap data
raw <- dataset.islet.rnaseq$raw[intersect.samples,]
samples <- dataset.islet.rnaseq$samples[intersect.samples,]
samples$DOwave <- as.factor(samples$DOwave)
samples$sex <- as.factor(samples$sex)



### Prefix to save name of .rds data
prefix <- 'attie_islet_mrna_284'

raw_file <- paste0(prefix,"_filtered_raw.rds")
norm_file <- paste0(prefix,"_normalized.rds")
norm_rz_file <- paste0(prefix,"_rZ_normalized.rds")
samples_file <-  paste0(prefix, "_samples_annot.rds")



### Normalize raw data by upper quantile for each subject
data.norm <- raw
q <- apply(data.norm, 1, quantile, probs=0.75)
data.norm <- data.norm/q                                  
data.norm <- log(data.norm + 1)



### Rank Z transformation function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


### Rankz across genes
data.rz <- data.norm
data.rz <- apply(data.rz, 2, rankZ)



### Saving the data to current working directory
saveRDS(raw, raw_file)
saveRDS(data.norm, norm_file)
saveRDS(data.rz, norm_rz_file)
saveRDS(samples, samples_file)





