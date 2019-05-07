library(clusterProfiler)
library(ggplot2)
library(ggraph)





load("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/Version 2/attie_islet_protein_hotspot_clusterprofiler.RData")






cp_data <- 'chr_2_164_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'nicely', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 2 Hotspot GO - Cellular Component') + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  geom_node_text(data = data, aes(x = x, y = y + .05), label = data$name, size = c(7,7,7,5,5,5,5,5,5))











cp_data <- 'chr_4_14_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'randomly', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Biological Processes') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
  geom_node_text(data = data, aes(x = x, y = y), label = data$name, size = c(7,7,5,5,5,5,5))




cp_data <- 'chr_4_14_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'sugiyama', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Cellular Component') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  geom_node_text(data = data, aes(x = x + c(.3,.3,.1,rep(.15,11)), y = y - c(-.02,.05,.24,rep(0,11))), label = data$name,
                 size = c(5,5,5,rep(5,11))) + 
  coord_flip()



cp_data <- 'chr_4_14_MF'
cp <- get(cp_data)
cp@result <- cp@result[-4,]
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot GO - Molecular Function') + theme(plot.title = element_text(hjust = 0.5, size = 20))






cp_data <- 'chr_4_14_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'kk', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 4 14 Hotspot KEGG Pathway') + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  geom_node_text(data  = data, aes(x = x, y = y), 
                 label = data$name,
                 size  = c(7,5,5,5,5,5))




















cp_data <- 'chr_5_139_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 139 Hotspot GO - Biological Processes') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
  geom_node_text(data  = data, aes(x = x, y = y), 
                 label = data$name,
                 size  = c(rep(5,8), rep(4,26)), repel = TRUE)




cp_data <- 'chr_5_139_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 139 Hotspot GO - Cellular Component') +
  theme(plot.title = element_text(hjust = 0.5, size = 20))




















cp_data <- 'chr_5_146_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'kk', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 5 146 Hotspot GO - Cellular Component') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
  geom_node_text(data  = data, aes(x = x, y = y), 
                 label = data$name,
                 size  = c(rep(7,2), rep(5,4)), repel = TRUE)



















cp_data <- 'chr_7_45_BP'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 82 Hotspot GO - Biological Processes') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))




cp_data <- 'chr_7_45_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 45 Hotspot GO - Cellular Component') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))




cp_data <- 'chr_7_45_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'kk', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 45 Hotspot GO - Molecular Function') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = 'none')






cp_data <- 'chr_7_45_KEGG'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'kk', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 45 Hotspot KEGG Pathway') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = 'none')






















cp_data <- 'chr_7_82_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'nicely', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + geom_node_text(aes_(x = data$x, y = data$y + c(rep(.03, 14), rep(0, nrow(data) - 14))), label = data$name, repel = TRUE) + theme(legend.position = 'none') +
  ggtitle('Chromsome 7 82 Hotspot GO - Cellular Component') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))




cp_data <- 'chr_7_82_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'nicely', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 7 82 Hotspot GO - Molecular Function') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = 'none')






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
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'nicely', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + theme(legend.position = 'none') +
  ggtitle('Chromsome 10 Hotspot GO - Biological Processes') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))



cp_data <- 'chr_10_41_CC'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = FALSE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + geom_node_text(aes_(x = data$x, y = data$y + c(rep(.03, 14), rep(0, nrow(data) - 14))), label = data$name, repel = TRUE) + theme(legend.position = 'none') +
  ggtitle('Chromsome 10 Hotspot GO - Cellular Component') +
  theme(plot.title = element_text(hjust = 0.5, size = 20))




cp_data <- 'chr_10_41_MF'
cp <- get(cp_data)
n <- nrow(cp)
plot_cp <- cnetplot(cp, showCategory = nrow(cp), layout = 'mds', colorEdge = TRUE, node_label = TRUE)
data <- plot_cp$data
data$color <- c(rep('cadetblue1', n), rep('darkseagreen1', nrow(plot_cp$data) - n))
plot_cp$data <- data
plot_cp + ggtitle('Chromsome 10 Hotspot GO - Molecular Function') +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = 'none')



