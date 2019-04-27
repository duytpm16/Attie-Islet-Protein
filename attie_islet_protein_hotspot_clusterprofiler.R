options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(plyr)
library(data.table)



load("~/Desktop/Attie Mass Spectrometry/QTL Viewers/attie_all_qtl_viewer_v5.RData")
lod.peaks <- dataset.islet.proteins$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
annots <- dataset.islet.proteins$annot.protein







## Counting number of cis LOD above 6 within a 4MB window across each chromosome
lod_df <- list()
slide <- 1
window <- 4
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(lod.peaks, qtl.chr == i)
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df     <- rbindlist(lod_df)






lod_df <- lod_df[lod_df$count > 45,]
lod_df <- lod_df %>% filter(chr %in% c('2','4','5','7','10','14') & pos %in% c(164, 14, 139, 146, 45, 82, 41, 101))












bg <- batchGenes(ids = unique(as.list(annots$gene.id)), version = 91)


for(i in 1:nrow(lod_df)){
  
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == lod_df$chr[i] & abs(qtl.pos - lod_df$pos[i]) <= 2)
  sub_lodpeaks$gene.id <- annots$gene.id[match(sub_lodpeaks$protein.id, annots$protein.id)]
  sub_lodpeaks$entrez.id <- bg$entrez_id[match(sub_lodpeaks$gene.id, bg$gene_id)]
  
  assign(paste0('chr_', lod_df$chr[i],'_',lod_df$pos[i],'_BP'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'BP', pAdjustMethod = 'fdr', readable = TRUE))

  assign(paste0('chr_', lod_df$chr[i],'_',lod_df$pos[i],'_CC'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'CC', pAdjustMethod = 'fdr', readable = TRUE))

  assign(paste0('chr_', lod_df$chr[i],'_',lod_df$pos[i],'_MF'),
         enrichGO(sub_lodpeaks$gene.id, universe = annots$gene.id, OrgDb = 'org.Mm.eg.db', keyType = 'ENSEMBL', ont = 'MF', pAdjustMethod = 'fdr', readable = TRUE))

  
  kegg <- enrichKEGG(sub_lodpeaks$entrez.id, universe = bg$entrez_id, organism = 'mmu', keyType = 'kegg', pAdjustMethod = 'fdr')
  for(j in 1:nrow(kegg@result)){
      kegg@result$geneID[j] <- paste0(bg$symbol[match(strsplit(kegg@result$geneID[j], split = '/')[[1]], bg$entrez_id)], collapse = '/')
  }
  assign(paste0('chr_', lod_df$chr[i],'_',lod_df$pos[i],'_KEGG'), kegg)
}







rm(list = ls()[!ls() %in% grep('chr_',ls(),value = TRUE)])
save.image("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Version 2/attie_islet_protein_hotspot_clusterprofiler.RData")

