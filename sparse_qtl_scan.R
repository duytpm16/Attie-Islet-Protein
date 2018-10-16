

rankz.protein <- dataset.islet.proteins$rankz
annot.protein <- dataset.islet.proteins$annots
covar.protein <- dataset.islet.proteins$covar




# Get chunk range
max_col = nrow(annot.protein)
rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
if(rng[length(rng)] > max_col) {
 rng = rng[1]:max_col
}


mediator.annot <- annot.protein
mediator.id <- mediator.annot$protein_id
mediator.chr <- mediator.annot$chr
mediator.start <- mediator.annot$start


target.annot <- annot.protein
target.id <- new.annot.protein$protein_id
target.chr <- new.annot.protein$chr
target.start <- new.annot.protein$start



sparse_data <- data.frame()
for(i in 1:nrow(annot.protein)){
  
  info_data <- as.data.frame(matrix(NA, nrow = nrow(mediator.annot) , ncol = 4))
  
  for(j in 1:nrow(annot.protein)){
      # Get overlapping samples in target and mediator data
      overlap.samples <- cbind(target = rankz.protein[,target.id[i], drop = FALSE], mediator = rankz.protein[,mediator.id[j], drop = FALSE])      
      overlap.samples <- overlap.samples[complete.cases(overlap.samples),]
      colnames(overlap.samples) <- c(target.id[i], mediator.id[j])
    
    
    
      # Find markers within 4Mbps of mediator start position
      temp <- subset(markers, markers$chr == mediator.chr[j] & abs(markers$pos - mediator.start[j]) <= 4)
    
    
    
      ### QTL scan on mediator markers
      gp = genoprobs[,mediator.chr[j]]
      gp[[1]] = gp[[1]][,,temp$marker, drop = FALSE]
      mediator.qtl <- scan1(gp, pheno = overlap.samples[,mediator.id[j], drop = FALSE], kinship = K[[mediator.chr[[j]]]], addcovar = covar.protein, cores = 0)
    
    
      # Find maximum marker position and ID
      max_qtl <- max_scan1(mediator.qtl, map = map)
    
    

      info_data[j,1:4] <- c(target.id[i], mediator.id[j], mediator.chr[j], rownames(max_qtl))
    
  }
  sparse_data <- bind_rows(sparse_data, info_data)
  print(i)
  
}

colnames(sparse_data) <- c('target.id', 'mediator.id', 'marker.id', 'target.lod', 'mediator.lod')

sparse_data[,c(4,5)] <- apply(sparse_data[,c(4,5)], 2, as.numeric)
saveRDS(sparse_data, file = paste0('attie_islet_protein_sparse_matrix_chunk_',args[1],'.rds'))

