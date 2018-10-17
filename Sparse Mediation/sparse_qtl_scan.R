library(dplyr)
library(intermediate2)
library(qtl2)

load('attie_islet_284_qtl_viewer.RData')





### Arguments. Number and size to break target and mediator into chunks
args = commandArgs(trailingOnly = TRUE)
targ.chunk_number = as.numeric(args[1])
targ.chunk_size = as.numeric(args[2])
med.chunk_number = as.numeric(args[3])
med.chunk_size = as.numeric(args[4])







### Extract rankZ, annotation, and covariate data
rankz.protein <- dataset.islet.proteins$rankz
annot.protein <- dataset.islet.proteins$annots
covar.protein <- dataset.islet.proteins$covar
stopifnot(colnames(rankz.protein) == annot.protein$protein_id)






### Chunk Range
# Get chunk range for mediator
max_col = nrow(annot.protein)
med.rng = ((med.chunk_number - 1) * med.chunk_size + 1):(med.chunk_number * med.chunk_size)
if(med.rng[length(med.rng)] > max_col) {
   med.rng = med.rng[1]:max_col
}

# Get chunk range for target
max_col = nrow(annot.protein)
targ.rng = ((targ.chunk_number - 1) * targ.chunk_size + 1):(targ.chunk_number * targ.chunk_size)
if(targ.rng[length(targ.rng)] > max_col) {
   targ.rng = targ.rng[1]:max_col
}







### Extract mediator and target data
# Mediator info
mediator.annot <- annot.protein[med.rng,]
mediator.id <- mediator.annot$protein_id
mediator.chr <- mediator.annot$chr
mediator.start <- mediator.annot$start

# Target info
target.annot <- annot.protein[targ.rng,]
target.id <- target.annot$protein_id
target.chr <- target.annot$chr
target.start <- target.annot$start







### Create empty data.frame to store target id, mediator id, mediator chromosome, mediator's best marker, 
#                                    sample size, target LOD and mediator LOD at mediator's best marker.
sparse_data <- as.data.frame(matrix(NA, nrow = nrow(target.annot) * nrow(mediator.annot) , ncol = 9))






### QTL scan begins
for(i in 1:nrow(target.annot)){
 
    # Create temporary index 
    index <- index <- ((i - 1) * med.chunk_size + 1):(i * med.chunk_size)
 
    for(j in 1:nrow(mediator.annot)){ 
        # Get overlapping samples in target and mediator data ( No NAs)
        overlap.samples <- cbind(target = rankz.protein[,target.id[i], drop = FALSE], 
                                 mediator = rankz.protein[,mediator.id[j], drop = FALSE])      
        overlap.samples <- overlap.samples[complete.cases(overlap.samples),]
        colnames(overlap.samples) <- c(target.id[i], mediator.id[j])
        sample.size <- nrow(overlap.samples)
    
    
        # Find markers within 4Mbps of mediator start position
        temp <- subset(markers, markers$chr == mediator.chr[j] & abs(markers$pos - mediator.start[j]) <= 4)
    
    
    
        ### QTL scan 
        # Mediator QTL scan
        gp = genoprobs[,mediator.chr[j]]
        gp[[1]] = gp[[1]][,,temp$marker, drop = FALSE]
        mediator.qtl <- scan1(gp, pheno = overlap.samples[,mediator.id[j], drop = FALSE], 
                              kinship = K[[mediator.chr[[j]]]], addcovar = covar.protein, cores = 4)
      
        # Find mediator's best marker (Highest LOD within 4Mbps of the start position)
        max_qtl <- max_scan1(mediator.qtl, map = map)
        max_marker <- rownames(max_qtl)
        gp[[1]] = gp[[1]][,,max_marker, drop = FALSE]
     
     
     
     
     
        # QTL scan on target at mediator's best marker
        target.qtl <- scan1(gp, pheno = overlap.samples[,target.id[i],drop = FALSE], kinship = K[[mediator.chr[j]]], addcovar = covar.protein, cores = 4)
    
        
     
     
     
     
        # Store LOD scores
        target.lod <- target.qtl[max_marker,]
        mediator.lod <- mediator.qtl[max_marker,]
    
    
    
     
     
     
        # Compute pearson and spearman correlation between the target and mediator
        pearson.cor <- cor(overlap.samples)[1,2]
        spearman.cor <- cor(overlap.samples, method = 'spearman')[1,2]
    
    
     
     
     
        # Store the results
        sparse_data[index[j],1:9] <- c(target.id[i], mediator.id[j], mediator.chr[j], max_marker, sample.size,
                                       target.lod, mediator.lod, pearson.cor, spearman.cor)
    
  }
  
  print(i)
  
}





### Fixing colnames and columns that should be numeric
colnames(sparse_data) <- c('target.id', 'mediator.id', 'chr', 'marker.id', 'sample.size',
                           'target.lod', 'mediator.lod', 'pearson', 'spearman')

sparse_data[,c(5,6,7,8,9)] <- apply(sparse_data[,c(5,6,7,8,9)], 2, as.numeric)


### Save data as .rds file
saveRDS(sparse_data, file = paste0('attie_islet_protein_sparse_matrix_chunk_',args[1],'_',args[3],'.rds'))

