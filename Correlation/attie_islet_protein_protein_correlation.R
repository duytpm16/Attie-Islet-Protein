### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)







### Load and Get Required Data
load("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Version 2/attie_islet_proteins_qtl_viewer_v2.RData")
full_protExpr   <- dataset.islet.proteins$rankz
full_protAnnots <- dataset.islet.proteins$annots










### Protein Protein correlation
prot_correlation <- cor(full_protExpr, full_protExpr, use = 'pairwise.complete.obs')











### Create dataframe to store correlation information
correlation_rownames <- rownames(prot_correlation)
correlation_colnames <- colnames(prot_correlation)
prot_df <- data.frame(prot.x.id     = correlation_colnames  [col(prot_correlation)],                                                               # Get protein id of protein x
                      prot.x.symbol = full_protAnnots$symbol[match(correlation_colnames, full_protAnnots$protein_id)][col(prot_correlation)],    # Get protein symbol of protein x
                      prot.x.chr    = full_protAnnots$chr   [match(correlation_colnames, full_protAnnots$protein_id)][col(prot_correlation)],       # Get protein chromosome of protein x
                      prot.x.start  = full_protAnnots$start [match(correlation_colnames, full_protAnnots$protein_id)][col(prot_correlation)],     # Get protein start of protein x
                      prot.x.end    = full_protAnnots$end   [match(correlation_colnames, full_protAnnots$protein_id)][col(prot_correlation)],       # Get protein end of protein x
                      prot.y.id     = rownames(prot_correlation)[row(prot_correlation)],                                                         # Get protein id of protein y
                      prot.y.symbol = full_protAnnots$symbol[match(correlation_rownames, full_protAnnots$protein_id)][row(prot_correlation)],    # Get protein symbol of protein y
                      prot.y.chr    = full_protAnnots$chr   [match(correlation_rownames, full_protAnnots$protein_id)][row(prot_correlation)],       # Get protein chromosome of protein y
                      prot.y.start  = full_protAnnots$start [match(correlation_rownames, full_protAnnots$protein_id)][row(prot_correlation)],     # Get protein start of protein y
                      prot.y.end    = full_protAnnots$end   [match(correlation_rownames, full_protAnnots$protein_id)][row(prot_correlation)],       # Get protein end of protein y
                      correlation   = c(prot_correlation))











### Arrange protein dataframe
prot_df <- prot_df %>% 
                    group_by(prot.x.id) %>%
                    arrange(prot.x.symbol, desc(correlation)) %>%
                    as.data.frame()











### Save as .rds
saveRDS(prot_df, '~/Desktop/attie_islet_protein_by_protein_correlation.rds')
