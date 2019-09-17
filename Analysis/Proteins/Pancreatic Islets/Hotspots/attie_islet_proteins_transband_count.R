### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)










### Load data
load('~/Desktop/Attie Mass Spectrometry/QTL Viewers/attie_all_qtl_viewer_v6.RData')
dataset <- 'dataset.islet.proteins'
peaks   <- get(dataset)$lod.peaks$additive











### Count number of QTLs within slide and window paramters. 
#     Note*: Warning is due to 'X' chromosome conversion to numeric. Need this to reorder 'chr' column.
transband_count <- function(markers, peaks, slide, window){

    counts_df <- markers %>% 
                  group_by(chr) %>% 
                  summarise(min = plyr::round_any(min(pos), accuracy = 1, f = floor),
                            max = plyr::round_any(max(pos), accuracy = window, f = ceiling),
                            len = length(seq(min, max, slide))) %>%
                  uncount(len) %>%
                  group_by(chr) %>%
                  mutate(pos = seq(unique(min), unique(max), slide), pos = ((pos + window) + pos)/2) %>%
                  select(-min, -max) %>%
                  group_by(chr, pos) %>%
                  mutate(counts = sum(peaks$qtl.chr == chr & abs(peaks$qtl.pos - pos) <= window /2)) %>%
                  arrange(as.numeric(chr), chr)
    counts_df             
}      








### Count number of qtls within tranband using LOD 6
peaks <- peaks %>% filter(!cis & lod > 6)
counts_lod6   <- transband_count(marker = markers, peaks = peaks, slide = 1, window = 4)


### Count number of qtls within tranband using LOD 7.3
peaks <- peaks %>% filter(!cis & lod > 7.3)
counts_lod7.3 <- transband_count(marker = markers, peaks = peaks, slide = 1, window = 4)










### Save
saveRDS(counts_lod6,   file = 'attie_islet_protein_transband_counts_lod6.rds')
saveRDS(counts_lod7.3, file = 'attie_islet_protein_transband_counts_lod7.3.rds')