library(dplyr)

### Load in data
#
#   1.) QTL_Viewer_V2.RData or Attie_DO378_eQTL_viewer_v6.Rdata
#   2.) qtl2 input data generated from qtl2_gather_mrna_scan1_input_data.R
#   3.) lodpeaks data: 37609 x 13
load("~/Desktop/Attie_QTL_Viewer_V2.RData")
load("~/Desktop/Islet Mediation/mRNA/attie_islet_mrna_284_qtl2_input.RData")
lod.peaks <- readRDS("attie_islet_mrna_284_rZ_lodpeaks_6.rds")


### Remove datasets in load data (1)
rm(dataset.islet.modules,dataset.cecum.metabolites,dataset.clinical.phenotypes,
   dataset.cecum.lipids,dataset.islet.proteins,dataset.liver.lipids,dataset.liver.metabolites,
   dataset.plasma.lipids,dataset.plasma.metabolites,dataset.islet.rnaseq.hotspots,dataset.islet.hotspots)



### Making changes to column names
#   Merging additional information to lodpeaks dataframe
#   Rename the merged info
colnames(lod.peaks) <- c('annot.id','marker.id','lod','qtl.chr','qtl.pos',LETTERS[1:8])
lod.peaks <- merge(lod.peaks, dataset.islet.rnaseq$annots[,c('gene_id','symbol','chr','start','end')], by.x = 'annot.id', by.y = 'gene_id', all.x = TRUE, )
lod.peaks <- lod.peaks %>% rename(gene.chr = chr,
                                  gene.start = start,
                                  gene.end = end,
                                  gene.symbol = symbol)


### Add cis column to lodpeaks
lod.peaks <- lod.peaks %>% mutate(cis = (qtl.chr == gene.chr & abs(qtl.pos - gene.start) <= 4))


### Reorder lodpeaks column
lod.peaks <- lod.peaks %>% select(annot.id, marker.id, lod, qtl.chr, qtl.pos, gene.symbol, gene.chr, gene.start, gene.end,cis, A,B,C,D,E,F,G,H)


### Put data into QTL Viewer format
dataset.islet.mrna <- list(annots = dataset.islet.rnaseq$annots,
                           covar = covar,
                           covar.factors = covar.factors,
                           datatype = datatype,
                           display.name = display.name,
                           expr = rankz,
                           norm = norm,
                           raw = raw,
                           lod.peaks = list(additive = lod.peaks, sex_int = data.frame()),
                           samples = samples)


### Remove data in environment that is in dataset.islet.mrna
rm(covar,covar.factors,datatype,display.name,rankz,norm,raw,samples,dataset.islet.rnaseq,lod.peaks)
