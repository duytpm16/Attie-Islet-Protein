load("~/Desktop/Attie/Attie_rZ_QTL_Viewer_V2.RData")


intersect.samples <- intersect(rownames(dataset.islet.rnaseq$expr), rownames(dataset.islet.proteins$expr))


raw <- dataset.islet.rnaseq$raw[intersect.samples,]
samples <- dataset.islet.rnaseq$samples[intersect.samples,]
samples$DOwave <- as.factor(samples$DOwave)
samples$sex <- as.factor(samples$sex)



prefix <- 'attie_islet_mrna_284'

raw_file <- paste0(prefix,"_filtered_raw.rds")
norm_file <- paste0(prefix,"_normalized.rds")
norm_rz_file <- paste0(prefix,"_rZ_normalized.rds")
samples_file <-  paste0(prefix, "_samples_annot.rds")



data.norm <- raw
# Normalized read counts in each sample
q <- apply(data.norm, 1, quantile, probs=0.75)
data.norm <- data.norm/q                                  






### Rank Z of normalized data.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


data.rz <- data.norm
data.rz <- apply(data.rz, 2, rankZ)



### Saving the data to current working directory
saveRDS(raw, raw_file)
saveRDS(data.norm, norm_file)
saveRDS(data.rz, norm_rz_file)
saveRDS(samples, samples_file)





