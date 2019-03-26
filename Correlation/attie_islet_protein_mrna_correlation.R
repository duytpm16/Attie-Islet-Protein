### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)







### Load and Get Required Data
load("~/Desktop/Attie Mass Spectrometry/Islet 284 /Version 2/attie_islet_284_qtl_viewer_v2.RData")
prot_expr   <- dataset.islet.proteins.284$rankz
prot_annots <- dataset.islet.proteins.284$annots
mrna_expr   <- dataset.islet.mrna.284$rankz
mrna_annots <- dataset.islet.mrna.284$annots
stopifnot(rownames(prot_expr) == rownames(mrna_expr))








### Protein mRNA correlation
correlation_matrix <- cor(mrna_expr, prot_expr, use = 'pairwise.complete.obs')











### Create dataframe to store correlation information
correlation_rownames <- rownames(correlation_matrix)
correlation_colnames <- colnames(correlation_matrix)
cor_df <- data.frame(prot.id     = correlation_colnames[col(correlation_matrix)],                                                  # Get protein id of protein x
                     prot.symbol = prot_annots$symbol[match(correlation_colnames, prot_annots$protein_id)][col(correlation_matrix)], # Get protein symbol of protein x
                     prot.chr    = prot_annots$chr   [match(correlation_colnames, prot_annots$protein_id)][col(correlation_matrix)], # Get protein chromosome of protein x
                     prot.start  = prot_annots$start [match(correlation_colnames, prot_annots$protein_id)][col(correlation_matrix)], # Get protein start of protein x
                     prot.end    = prot_annots$end   [match(correlation_colnames, prot_annots$protein_id)][col(correlation_matrix)], # Get protein end of protein x
                     mrna.id     = rownames(correlation_matrix)[row(correlation_matrix)],                                            # Get protein id of protein y
                     mrna.symbol = mrna_annots$symbol[match(correlation_rownames, mrna_annots$gene_id)][row(correlation_matrix)],    # Get protein symbol of protein y
                     mrna.chr    = mrna_annots$chr   [match(correlation_rownames, mrna_annots$gene_id)][row(correlation_matrix)],    # Get protein chromosome of protein y
                     mrna.start  = mrna_annots$start [match(correlation_rownames, mrna_annots$gene_id)][row(correlation_matrix)],    # Get protein start of protein y
                     mrna.end    = mrna_annots$end   [match(correlation_rownames, mrna_annots$gene_id)][row(correlation_matrix)],    # Get protein end of protein y
                     correlation   = c(correlation_matrix))











### Arrange protein dataframe
cor_df <- cor_df %>% 
                   group_by(prot.id) %>%
                   arrange(prot.symbol, desc(correlation)) %>%
                   as.data.frame()






cor_df <- setkey(as.data.table(cor_df), prot.symbol, correlation)




### Save as .rds
saveRDS(cor_df, '~/Desktop/attie_islet_protein_by_protein_correlation.rds')
