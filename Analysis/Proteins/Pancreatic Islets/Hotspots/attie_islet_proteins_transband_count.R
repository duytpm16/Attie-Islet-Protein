### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)










### Load data
load('~/Desktop/Attie Mass Spectrometry/QTL Viewers/attie_all_qtl_viewer_v6.RData')
dataset <- 'dataset.islet.proteins'
peaks   <- get(dataset)$lod.peaks$additive










### Help function:
#       Counts number of QTLs within certain window across genome.
#       
#       Parameters:
#           markers - marker's dataframe
#           peaks   - QTL peaks dataframe with 'qtl.chr' and 'qtl.pos' at least in the column names
#           slide   - Mbp to slide across genome
#           window  - Width of window in Mbps
#
#       Note*: 
#           Warning is due to 'X' chromosome conversion to numeric. Need this to reorder 'chr' column (line 52).
transband_count <- function(markers, peaks, slide, window){
  
      markers %>% 
          group_by(chr) %>% 
          summarise(min = plyr::round_any(min(pos), accuracy = 1, f = floor),
                    max = plyr::round_any(max(pos), accuracy = window, f = ceiling),
                    len = length(seq(min, max, slide))) %>%
          uncount(len) %>% group_by(chr) %>%
          mutate(pos = seq(unique(min), unique(max), slide), pos = ((pos + window) + pos)/2) %>%
          select(-min, -max) %>% group_by(chr, pos) %>%
          mutate(counts = sum(peaks$qtl.chr == chr & abs(peaks$qtl.pos - pos) <= window /2)) %>%
          arrange(as.numeric(chr), chr)
  
}      













### Count number of qtls with LOD 6 in 4Mbp windows and 1 Mbp step across genome. Finally select windows with highest counts
peaks6 <- peaks %>% filter(!cis & lod > 6)
counts_lod6   <- transband_count(marker = markers, peaks = peaks6, slide = 1, window = 4)
counts_lod6   <- counts_lod6 %>% 
                    filter(counts >= 45) %>% 
                    filter(chr %in% c('2', '4', '5', '7', '10', '14')) %>%
                    filter(pos %in% c(164, 14, 139, 146, 45, 82, 41, 101))


### Count number of qtls with LOD 7.3 in 4Mbp windows and 1 Mbp step across genome. Finally select windows with highest counts
peaks7.3 <- peaks %>% filter(!cis & lod > 7.3)
counts_lod7.3 <- transband_count(marker = markers, peaks = peaks7.3, slide = 1, window = 4)
counts_lod7.3 <- counts_lod7.3 %>% 
                    filter(counts >= 15) %>% 
                    filter(chr %in% c('2', '4', '5', '7')) %>%
                    filter(pos %in% c(163, 13, 139, 146, 45, 82))













### Save
saveRDS(counts_lod6,   file = 'attie_islet_protein_transband_counts_lod6.rds')
saveRDS(counts_lod7.3, file = 'attie_islet_protein_transband_counts_lod7.3.rds')
