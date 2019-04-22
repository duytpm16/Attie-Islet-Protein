## Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)





### Load and get data
#load("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/attie_islet_protein_rZ_qtl_viewer.RData")
lod.peaks    <- dataset.islet.proteins$lod.peaks$additive
lod.peaks    <- lod.peaks[complete.cases(lod.peaks),]
lod.peaks    <- lod.peaks[lod.peaks$gene.chr %in% c(1:19,'X'),]
local.peaks  <- lod.peaks[lod.peaks$cis,]
distal.peaks <- lod.peaks[!lod.peaks$cis,]
colnames(lod.peaks) <- gsub('.','_',colnames(lod.peaks), fixed = TRUE)
dis_color <- 'deepskyblue4'
cis_color <- 'lightgreen'




### Format chromosomes for transcriptome plot
all.chr = lod.peaks %>%
  select(qtl_chr, gene_chr) %>%
  gather(k, v) %>%
  select(v) %>%
  distinct() %>%
  arrange(v)

all.chr <- all.chr$v[!is.na(all.chr$v)]

if(length(grep("M", all.chr)) > 0){
  wh <- grep("M", all.chr)
  all.chr <- all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
}




### Format lod.peaks for transcriptome plot
data = lod.peaks %>% 
  mutate(cis      = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= 4),
         qtl_chr  = factor(qtl_chr,  
                           levels   = all.chr[order(as.numeric(all.chr))]),
         gene_chr = factor(gene_chr,  
                           levels   = rev(all.chr[order(as.numeric(all.chr))])),
         gene_pos = (gene_end + gene_start) * 0.5)
cis.colors        = c(dis_color, cis_color)
names(cis.colors) = c("FALSE", "TRUE")










## Counting number of cis LOD above 6 within a 4MB window across each chromosome
lod_df_cis <- list()
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
  sub <- subset(local.peaks, qtl.chr == i)
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_cis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}
rm(min,max,x_axis_scale,chr,sub,count,pos)



lod_df_cis     <- rbindlist(lod_df_cis)
lod_df_cis$chr <- factor(lod_df_cis$chr, levels = c(1:19,'X'))










## Counting number of distal LOD above 6 within a 4MB window across each chromosome
lod_df_dis <- list()
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
  sub <- subset(distal.peaks, qtl.chr == i)
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_dis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}
rm(min,max,x_axis_scale,chr,sub,count,pos)



lod_df_dis <- rbindlist(lod_df_dis)
lod_df_dis$chr <- factor(lod_df_dis$chr, levels = c(1:19,'X'))










### Plot cis density plot
c <- ggplot(lod_df_cis, aes(x = pos, y = count)) +
  geom_line(col = cis_color) +
  geom_hline(yintercept = 45, linetype = 2, color = 'grey70') +
  ylab('No. local pQTL/4 Mbp') +
  xlab('Chromosome') +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.margin = unit(c(0,0,1.3,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "grey70"),
        axis.text.x = element_blank(),
        strip.text = element_text(size= 23),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
  facet_grid(.~chr, scales = "free",switch = 'x')








### Plot distal density plot
d <- ggplot(lod_df_dis, aes(x = pos, y = count)) +
  geom_line(col = dis_color) +
  geom_hline(yintercept = 45, linetype = 2, color = 'grey70') +
  ylab('No. distal pQTL/4 Mbp') +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.margin = unit(c(2,0,0,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "grey70"),
        axis.text.x = element_blank(),
        strip.text = element_text(size= 23),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
  facet_grid(.~chr, scales = "free")













### Plot transcriptome
g <- ggplot(data, aes(x = qtl_pos, y = gene_pos, size = log(lod)), alpha = 0.4) +
  geom_point(aes(color = cis), alpha = 0.4) + 
  scale_color_manual(values = cis.colors, labels = c('Distal-pQTL','Local-pQTL')) +
  facet_grid(gene_chr ~ qtl_chr, scales = "free",switch = 'y', shrink = TRUE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "grey70"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.0, "line"),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 180),
        axis.text = element_blank(),
        axis.line.y.left = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_blank(),
        legend.key = element_blank()) +
  xlab('QTL position (Chr)') +
  ylab('Gene position (Chr)')








### Put plots in 3 panels vertically. DONT LOAD COWPLOT LIBRARY!!
cowplot::plot_grid(g, d,c, axis = 'lr', ncol = 1, align = 'v', rel_heights = c(1,.3,.3), labels = c('A','B','C'), vjust = c(1.5,7,1.5))

#grid.text(label = 'Distal-pQTL', x=unit(.15,'npc'), y = unit(.97,'npc'),gp=gpar(fontsize=20, fontface = 'bold', col = dis_color))
#grid.text(label = 'Local-pQTL', x=unit(.15,'npc'), y = unit(.95,'npc'), gp=gpar(fontsize=20, fontface = 'bold',col=cis_color))



### Put arrows on top of hotspots
# Chromosome 1
grid.lines(x = unit(c(0.09, .09), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 2
grid.lines(x = unit(c(0.135, .135), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 4
grid.lines(x = unit(c(0.195, .195), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 5
grid.lines(x = unit(c(0.277, .277), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))

grid.lines(x = unit(c(0.28, .28), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 7
grid.lines(x = unit(c(0.345, .345), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))
grid.lines(x = unit(c(0.355, .355), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 10
grid.lines(x = unit(c(0.487, .487), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))


# Chromosome 14
grid.lines(x = unit(c(0.698, .698), "npc"),
           y = unit(c(0.33, .329), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"),type="closed"))
