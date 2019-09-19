### Options and libraries
options(stringsAsFactors = FALSE)
library(VennDiagram)
library(tidyverse)
library(reshape2)
library(cowplot)






### Load  and get data
load('attie_all_qtl_viewer_v6.RData')
prot <- 'dataset.islet.proteins'
mrna <- 'dataset.islet.rnaseq'


prot_annots <- get(prot)[['annot.protein']]
prot_peaks  <- get(prot)[['lod.peaks']]$additive %>% 
                  merge(., prot_annots[,c('protein.id','gene.id')], all.x = TRUE, sort = FALSE) %>% 
                  filter(lod > 7.2)

mrna_annots <- get(mrna)[['annot.mrna']]
mrna_peaks  <- get(mrna)[['lod.peaks']]$additive %>% filter(lod > 7.2)












### Find proteins with a transcript measurement
prot_annots <- prot_annots %>% filter(gene.id %in% mrna_annots$gene.id)
mrna_annots <- mrna_annots %>% filter(gene.id %in% prot_annots$gene.id)











### Find number of local pQTLs with local eQTLs
prot_local <- prot_peaks %>% filter(cis) %>% filter(gene.id %in% prot_annots$gene.id)
mrna_local <- mrna_peaks %>% filter(cis) %>% filter(gene.id %in% mrna_annots$gene.id)
n_local <- sum(prot_local$gene.id %in% mrna_local$gene.id)









### Find number of distal pQTLs with distal eQTLs at the same location (within 4 Mbp)
prot_distal <- prot_peaks %>% filter(!cis) %>% filter(gene.id %in% prot_annots$gene.id)
mrna_distal <- mrna_peaks %>% filter(!cis) %>% filter(gene.id %in% mrna_annots$gene.id)

distal_overlaps <- prot_distal %>%
                      group_by(protein.id, gene.id, gene.symbol, qtl.chr, qtl.pos) %>%
                      summarise(count = sum(mrna_distal$qtl.chr == qtl.chr & 
                                            abs(mrna_distal$qtl.pos - qtl.pos) <= 4 &
                                            mrna_distal$gene.id == gene.id))
n_distal <- sum(distal_overlaps$count)










### Create table of QTL overlap info
result <- data.frame(protein = c(nrow(prot_local), nrow(prot_distal)),
                     mrna    = c(nrow(mrna_local), nrow(mrna_distal)),
                     overlap = c(n_local, n_distal))
rownames(result) <- c('local', 'distal')














## local venn diagram
local_venn <- draw.pairwise.venn(area1      = result['local','protein'],
                                 area2      = result['local','mrna'],
                                 cross.area = result['local','overlap'],
                                 scaled     = TRUE,
                                 category   = c('Proteins', 'Transcripts'),
                                 cat.pos    = c(180, 120),
                                 cat.just   = list(c(.6,7.5), c(-.55,1.2)),
                                 cat.cex    = 1.5,
                                 cat.fontface = 2,
                                 col        = c('dodgerblue3','indianred3'),
                                 label.col  = c('grey30'),
                                 fontface   = 1,
                                 lwd        = 4,
                                 lty        = 6,
                                 ext.line.lwd = 3,
                                 ext.line.lty = 2,
                                 ext.pos    = 90,
                                 ext.length = .5,
                                 fill       = c("grey97"),
                                 cex        = c(2.,2.,2.5),
                                 rotation.degree = -45)





## distal venn diagram
distal_venn <- draw.pairwise.venn(area1      = result['distal','protein'],
                                  area2      = result['distal','mrna'],
                                  cross.area = result['distal','overlap'],
                                  scaled     = TRUE,
                                  category   = c('Proteins', 'Transcripts'),
                                  cat.pos    = c(180, 120),
                                  cat.just   = list(c(.45,9.5), c(-.55,3)),
                                  cat.cex    = 1.5,
                                  cat.fontface = 2,
                                  col        = c('dodgerblue3','indianred3'),
                                  label.col  = c('grey30'),
                                  fontface   = 1,
                                  lwd        = 4,
                                  lty        = 6,
                                  ext.line.lty = 2,
                                  ext.line.lwd = 3,
                                  ext.pos    = 0,
                                  ext.length = .5,
                                  fill       = c("grey97"),
                                  cex        = c(2.,2.,2.5),
                                  rotation.degree = -45)
                     
            









### Plot both vennDiagram together
plot_grid(grobTree(local_venn), grobTree(distal_venn),
          ncol = 2, 
          labels = c('Local QTLs','Distal QTLs'), 
          label_size = 18, label_fontface = 2, rel_heights = c(.1,.8,.8), 
          hjust = c(0,.5))


