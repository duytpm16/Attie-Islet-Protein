############################################################################################################################
#
#   This script takes the LOD peaks dataframe to find hotspots and plot a heatmap of the allele effect of each
#       expressions within the hotspot.
#
#   Note*: Will have to manually select the hotspot
#
#
#   Input:
#      1.) QTL Viewer dataset
#      2.) Which dataset to grab
#      3.) Window of hotspot
#      4.) Slide in MB across genome
#
#
#   Output:
#      1.) Heatmap plot of allele effects for each expression in the hotspot
#
#
#  
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#   Date:   November 26, 2018
#############################################################################################################################
options(stringsAsFactors = F)
options(scipen = 999)


library(plyr)
library(data.table)
library(tidyverse)
library(gplots)


### Variables to change
load("~/Desktop/Attie Not Imputed and Sex+Wave/DO 284/attie_islet_284_nobatch_qtl_viewer.RData")
dataset   <- 'dataset.islet.proteins'
window    <- 4
slide     <- 1







### Extract required data
pheno            <- get(dataset)$rankz
lod.peaks        <- get(dataset)$lod.peaks$additive












### Counting number of LOD above 6 within a 4MB window across each chromosome
lod.peaks <- lod.peaks[lod.peaks$cis == FALSE,]
lod_df <- list()
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


lod_df <- rbindlist(lod_df)











### Select chromosomes with large hotspot peak
hotspot <- lod_df[which(lod_df$count >= quantile(lod_df$count, .98)[[1]]),]














### Color for heatmap
hmcols<-colorRampPalette(c("blue","white","red"))(256)





### Heatmap plot begin
#    1.) Get proteins within selected hotspots
#    2.) Order rows by expression correlation
#    3.) Plot heatmap of allele effects
for(i in c(3,8,10,14,22,25,34,38,46,51,54)){

    # 1
    sub <- subset(lod.peaks, qtl.chr == hotspot$chr[i] & abs(qtl.pos - hotspot$pos[i]) <= window/2 & cis == FALSE)

    
    
    # 2
    expr.cor = cor(pheno[,sub$annot.id], use = 'pairwise.complete.obs')
    cl = hclust(as.dist(1.0 - expr.cor), method = "average")
    sub = sub[cl$order,]
    rownames(sub) <- sub$gene.symbol

    
    
    # 3
    png(paste0('attie_islet_proteins_',hotspot$chr[i],'_',hotspot$pos[i],'_allele_effect_heatmap.png'), width = 960, 960)
    heatmap.2(as.matrix(sub[,LETTERS[1:8]]), Rowv = FALSE, Colv = FALSE, dendrogram = 'none',
              trace = 'none', colsep = 1:8, sepwidth = c(0,.1), sepcolor = 'black',
              col = hmcols, offsetCol = -88, offsetRow = -0.25,
              srtCol = 360, margins = c(8,10), cexRow = 1,
              cellnote = round(as.matrix(sub[,LETTERS[1:8]]),4), notecol = 'black', notecex = 1.2,
              cexCol = 2, main = paste('\n\n\nChromosome', hotspot$chr[i],'  Position', hotspot$pos[i]))
    dev.off()

}
  

