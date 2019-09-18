### Options and libraries
options(stringsAsFactors = FALSE)
library(clusterProfiler)
library(tidyverse)







### Load data
load('attie_all_qtl_viewer_v6.RData')
dataset <- 'dataset.islet.proteins'
annots  <- get(dataset)[['annot.protein']]
peaks   <- get(dataset)[['lod.peaks']]$additive
transband6   <- readRDS('attie_islet_protein_transband_counts_lod6.rds')
transband7.3 <- readRDS('attie_islet_protein_transband_counts_lod7.3.rds')







### Remove cis and add Ensembl Gene IDs to peaks
peaks <- peaks %>% filter(!cis)
peaks <- merge(x  = peaks, y = annots[,c('protein.id', 'gene.id')],
               by = 'protein.id', all.x = TRUE, sort = FALSE)










### Help function:
#     Performs enrichment for each transband on ENSEMBL ids. 
#       Using a distance of +/- 2 from pos since I use 4Mbp to define transband.
#
#     Parameters:
#       transband_df - dataframe containing 'chr', 'pos', and 'counts' columns
#       peaks_df     - peaks dataframe with at least 'qtl.chr' and 'qtl.pos' columns
#       annots_df    - annotation dataframe containing all ids for every expression.
#                         Needs at least 'gene.id' column
transband_enrichment <- function(transband_df, peaks_df, annots_df){
  
    prot_list <- list()
    enrich_bp <- list()
    enrich_cc <- list()
    enrich_mf <- list()
    for(i in 1:nrow(transband_df)){
      
        # Defining name for list indices
        name <- paste0('chr',transband_df$chr[i], '_', transband_df$pos[i])
        
        
        # Subset peaks in transband
        subPeak <- peaks_df %>%
                     filter(qtl.chr == transband_df$chr[i] & 
                            abs(qtl.pos - transband_df$pos[i]) <= 2)
        stopifnot(nrow(subPeak) == transband_df$counts[i])
    
        
        
        
        
        
        # Enrichment on proteins in transband
        bp <- enrichGO(gene = subPeak$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL',
                       ont  = 'BP', universe = unique(annots_df$gene.id), pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
        
        cc <- enrichGO(gene = subPeak$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL',
                       ont  = 'CC', universe = unique(annots_df$gene.id), pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
        
        mf <- enrichGO(gene = subPeak$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL',
                       ont  = 'MF', universe = unique(annots_df$gene.id), pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
        
        
        
        
        
        
        # Store results
        prot_list[[name]] <- subPeak
        enrich_bp[[name]] <- bp
        enrich_cc[[name]] <- cc
        enrich_mf[[name]] <- mf
        
        
        
        print(paste0(i, ' of ', nrow(transband_df)))
    }
    
    
    
    
    # Return list
    return(list('prot_list' = prot_list, 
                'enrich_bp' = enrich_bp, 
                'enrich_cc' = enrich_cc, 
                'enrich_mf' = enrich_mf))
}











### Run enrichment function for different LOD thresholds
peaks6 <- peaks %>% filter(lod > 6)
results_lod6   <- transband_enrichment(transband_df = transband6, peaks_df = peaks6, annots_df = annots)

peaks7.3 <- peaks %>% filter(lod > 7.3)
results_lod7.3 <- transband_enrichment(transband_df = transband7.3, peaks_df = peaks7.3, annots_df = annots)













### Save 
prot_list <- results_lod6[['prot_list']]
enrich_bp <- results_lod6[['enrich_bp']]
enrich_cc <- results_lod6[['enrich_cc']]
enrich_mf <- results_lod6[['enrich_mf']]
save(prot_list, enrich_bp, enrich_cc, enrich_mf, file = 'attie_islet_proteins_transband_enrichment_lod6.Rdata')




prot_list <- results_lod7.3[['prot_list']]
enrich_bp <- results_lod7.3[['enrich_bp']]
enrich_cc <- results_lod7.3[['enrich_cc']]
enrich_mf <- results_lod7.3[['enrich_mf']]
save(prot_list, enrich_bp, enrich_cc, enrich_mf, file = 'attie_islet_proteins_transband_enrichment_lod7.3.Rdata')



