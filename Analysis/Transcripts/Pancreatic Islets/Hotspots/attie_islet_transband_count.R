### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)










### Load data
load('~/Desktop/Attie Mass Spectrometry/QTL Viewers/attie_all_qtl_viewer_v6.RData')
dataset <- 'dataset.islet.rnaseq'
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















### Count number of qtls with LOD 7.3 in 4Mbp windows and 1 Mbp step across genome
peaks7.2 <- peaks %>% filter(!cis & lod > 7.2)
counts_lod7.2 <- transband_count(marker = markers, peaks = peaks7.2, slide = 1, window = 4)
counts_lod7.2 <- counts_lod7.2 %>% 
                    filter(counts >= 100) %>%
                    filter(chr %in% c('2','5','7','11','13')) %>%
                    filter(pos %in% c(164, 146, 46, 71, 111))













### Save
saveRDS(counts_lod7.2, file = 'attie_islet_transcript_transband_counts_lod7.2.rds')
