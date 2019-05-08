options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(intermediate)





load('~/Desktop/Attie Mass Spectrometry/QTL Viewers/attie_all_qtl_viewer_v5.RData')
target_ds   <- 'dataset.islet.proteins'
mediator_ds <- 'dataset.islet.rnaseq'
target_id   <- 'protein.id'
mediator_id <- 'gene.id'

overlap <- intersect(get(target_ds)$annot.samples$mouse.id, get(mediator_ds)$annot.samples$mouse.id)

lod.peaks <- get(target_ds)$lod.peaks$additive
lod.peaks <- lod.peaks[!lod.peaks$cis,]
expr      <- get(target_ds)$data$rz[overlap,]
covar     <- get(target_ds)$covar.matrix[overlap,,drop = FALSE]



annot <- switch(get(mediator_ds)$datatype,"mRNA" = 'annot.mrna', "mrna" = 'annot.mrna', "protein" = 'annot.protein') 
med_expr  <- get(mediator_ds)$data$rz[overlap,]
med_annot <- get(mediator_ds)[[annot]] %>% dplyr::rename(pos = start)
med_expr  <- med_expr[rownames(expr), med_annot[,mediator_id, drop = TRUE]]
med_type  <- get(mediator_ds)$datatype






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
lod_df$chr <- factor(lod_df$chr, levels = c(1:19,'X'))





lod_df <- lod_df[lod_df$count > 45,]
lod_df <- lod_df %>% filter(chr %in% c('2','4','5','7','10','14') & pos %in% c(164, 14, 139, 146, 45, 82, 41, 101))














interval <- 2
color <- 'dodgerblue1'
for(i in 1:nrow(lod_df)){
  
  hs_chr <- as.character(lod_df$chr[i])
  
  sub_lodpeaks <- lod.peaks %>% subset(qtl.chr == lod_df$chr[i] & abs(qtl.pos - lod_df$pos[i]) <= 2)
  sub_expr     <- expr[,sub_lodpeaks[,target_id, drop = TRUE]]
  stopifnot(sub_lodpeaks[,target_id] == colnames(sub_expr))
  
  
  med_in_interval <- med_annot %>% subset(chr == lod_df$chr[i] & abs(pos - lod_df$pos[i]) <= interval) %>% dplyr::arrange(pos)
  med_in_interval_expr <- med_expr[,med_in_interval[,mediator_id, drop = TRUE]]
  stopifnot(colnames(med_in_interval_expr) == med_in_interval[,mediator_id, drop = FALSE])
  
  
  
  
  corr  <- cor(sub_expr, use = 'pairwise.complete.obs')
  cl    <- hclust(as.dist(1.0 - corr), method = "average") 
  order <- sub_lodpeaks[cl$order,]
  med_result  <- as.data.frame(matrix(0, nrow = nrow(order), ncol = ncol(med_in_interval_expr),dimnames = list(order[,target_id, drop = TRUE], colnames(med_in_interval_expr))))
  stopifnot(rownames(sub_expr) == rownames(med_in_interval_expr))
  stopifnot(colnames(med_in_interval_expr) == colnames(med_result))
  stopifnot(rownames(sub_expr) == rownames(covar))
  
  
  
  
  
  
  for(j in 1:nrow(order)){
    med <- mediation.scan(target     = sub_expr[, order[,target_id, drop = TRUE][j], drop = FALSE],
                          mediator   = med_in_interval_expr,
                          qtl.geno   = genoprobs[[order$qtl.chr[j]]][rownames(sub_expr),,order$marker.id[j]],
                          annotation = med_in_interval,
                          covar      = covar,
                          method     = 'double-lod-diff',
                          verbose    = FALSE)
    
    med$dp <- (order$lod[j] - med$LOD) / order$lod[j]
    med$dp[med$dp < 0] <- 0
    med_result[order[,target_id, drop = TRUE][j], med[,mediator_id]] <- med$dp
    
  }
  stopifnot(rownames(med_result) == order[,target_id, drop = TRUE])
  stopifnot(colnames(med_result) == med_in_interval[,mediator_id, drop = TRUE])
  
  
  
  
  
  
  order$gene.symbol[duplicated(order$gene.symbol)] <- paste0(order$gene.symbol[duplicated(order$gene.symbol)],'-1')
  med_in_interval$symbol[duplicated(med_in_interval$symbol)] <- paste0(med_in_interval$symbol[duplicated(med_in_interval$symbol)],'-1')
  rownames(med_result) <- order$gene.symbol
  colnames(med_result) <- med_in_interval$symbol
  
  
  
  
  
  med_melt <- cbind(target = rownames(med_result), med_result)
  med_melt <- melt(med_melt)
  med_melt$target <- factor(med_melt$target, levels = rownames(med_result))
  med_melt$variable <- factor(med_melt$variable, levels = colnames(med_result))
  med_melt <- med_melt %>% plyr::rename(c("value" = "Proportion Drop"))
  
  
  
  
  
  
  
  print(ggplot(med_melt, aes(x = variable, y = target)) + geom_tile(aes(fill = med_melt$`Proportion Drop`), color = 'black')  + 
          xlab(paste0('Chr',hs_chr, '  ~ ', lod_df$pos[i] - interval, '-', lod_df$pos[i] + interval, ' Mbp')) + 
          ylab(paste0('Proteins in Chr',hs_chr,' HotSpot')) +
          theme(plot.title = element_text(hjust = .5, size = 18),
                axis.title = element_text(size= 18),
                axis.text.x = element_text(angle = 90, vjust = .5, face = 'bold', size = 12),
                axis.text.y = element_text(vjust = .5, face = 'bold', size = 12),
                panel.background = element_blank(),
                legend.title = element_blank()) +
          scale_fill_gradientn(colors = c(color,'orange','yellow')) + ggtitle(paste0('Chr',hs_chr,' HotSpot Interval Mediation')))
  
  
  
}

