library(clusterProfiler)
library(ggplot2)






load("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Version 2/attie_islet_protein_hotspot_clusterprofiler.RData")






cp_data <- 'chr_2_164_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 2 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))











cp_data <- 'chr_4_14_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_4_14_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_4_14_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_4_14_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5))




















cp_data <- 'chr_5_139_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 139 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_5_139_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 139 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




















cp_data <- 'chr_5_146_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 146 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))






















cp_data <- 'chr_7_45_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 82 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_7_45_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 45 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_7_45_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 45 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')






cp_data <- 'chr_7_45_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 45 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')






















cp_data <- 'chr_7_82_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + geom_node_text(aes_(x = data$x, y = data$y + c(rep(.03, 14), rep(0, nrow(data) - 14))), label = data$name, repel = TRUE) + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 82 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_7_82_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 82 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')






cp_data <- 'chr_7_82_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 82 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')











cp_data <- 'chr_10_41_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 10 Hotspot GO - Biological Processes') + theme(plot.title = element_text(hjust = 0.5))



cp_data <- 'chr_10_41_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + geom_node_text(aes_(x = data$x, y = data$y + c(rep(.03, 14), rep(0, nrow(data) - 14))), label = data$name, repel = TRUE) + theme(legend.position = 'none') +
  ggtitle('Chromsome 10 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5))




cp_data <- 'chr_10_41_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 10 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')





