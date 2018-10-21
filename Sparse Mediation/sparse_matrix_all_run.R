library(dplyr)
library(tidyr)
library(intermediate2)
library(qtl2)

load('attie_islet_284_qtl_viewer.RData')





### Arguments. Number and size to break target and mediator into chunks
args = commandArgs(trailingOnly = TRUE)
targ.chunk_number = as.numeric(args[1])
targ.chunk_size = as.numeric(args[2])
med.chunk_number = as.numeric(args[3])
med.chunk_size = as.numeric(args[4])
targ = args[5]
med = args[6]
targ_id = args[7]
med_id = args[8]





### Extract rankZ, annotation, and covariate data
rankz.targ <- get(targ)$rankz
rankz.med <- get(med)$rankz
annot.targ <- get(targ)$annots
annot.med <- get(med)$annots
annot.targ <- annot.targ %>% dplyr::rename(id = targ_id)
annot.med <- annot.med %>% dplyr::rename(id = med_id)
annot.med$chr <- as.character(annot.med$chr)
annot.med <- annot.med %>% filter(chr %in% c(1:19,'X'))
rownames(annot.med) <- annot.med$id
rankz.med <- rankz.med[,annot.med$id]
covar.targ <- get(targ)$covar
covar.med <- get(med)$covar
stopifnot(colnames(rankz.targ) == annot.targ$id)
stopifnot(colnames(rankz.med) == annot.med$id)







### Chunk Range
# Get chunk range for mediator
max_col = nrow(annot.med)
med.rng = ((med.chunk_number - 1) * med.chunk_size + 1):(med.chunk_number * med.chunk_size)
if(med.rng[length(med.rng)] > max_col) {
  med.rng = med.rng[1]:max_col
}

# Get chunk range for target
max_col = nrow(annot.targ)
targ.rng = ((targ.chunk_number - 1) * targ.chunk_size + 1):(targ.chunk_number * targ.chunk_size)
if(targ.rng[length(targ.rng)] > max_col) {
  targ.rng = targ.rng[1]:max_col
}







### Extract mediator and target data
# Mediator info
mediator.annot <- annot.med[med.rng,]

# Target info
target.annot <- annot.targ[targ.rng,]




### Create vectors to store the data

#     Create dataframe containing all possible combinations of target and mediator. Then add mediator's chromosome and start position
target_and_mediator = crossing(target.annot$id, mediator.annot$id)
target_and_mediator$chr = mediator.annot$chr[match(target_and_mediator$`mediator.annot$id`, mediator.annot$id)]
target_and_mediator$start = mediator.annot$start[match(target_and_mediator$`mediator.annot$id`, mediator.annot$id)]


#     Create vector of the ids 
target.id <- target_and_mediator$`target.annot$id`
mediator.id <- target_and_mediator$`mediator.annot$id`
mediator.chr <- target_and_mediator$chr
mediator.start <- target_and_mediator$start


#     Initialize empty vectors to store info
n <- nrow(target_and_mediator)
best.marker <- character(length = n)
sample.size <- numeric(length = n)
target.lod <- numeric(length = n)
mediator.lod <- numeric(length = n)
pearson <- numeric(length = n)
spearman <- numeric(length = n)
best.model <- character(length = n)
best.model.p <- numeric(length = n)
best.triad <- character(length = n)
best.triad.p <- numeric(length = n)
mediation.lod <- numeric(length = n)
inverse.mediation.lod <- numeric(length = n)





### QTL scan and med begins
for(i in 1:nrow(target_and_mediator)){
  
    # Create temporary index 
    target.name <- target.id[i]
    current.mediator <- mediator.id[i]
    current.mediator.chr <- mediator.chr[i]
    kin <- K[[current.mediator.chr]]
      
      
    # Get overlapping samples in target and mediator data ( No NAs)
    overlap.samples <- cbind(rankz.targ[,target.name,drop = FALSE], rankz.med[,current.mediator,drop = FALSE])      
    overlap.samples <- overlap.samples[complete.cases(overlap.samples),]
      
      
    # Find markers within 4Mbps of mediator start position
    temp <- subset(markers, markers$chr == current.mediator.chr & abs(markers$pos - mediator.start[i]) <= 4)
      
    ### QTL scan on mediator markers
    gp = genoprobs[,current.mediator.chr] 
    gp[[1]] = gp[[1]][,,temp$marker, drop = FALSE]
    mediator.qtl <- scan1(gp, pheno = overlap.samples[,current.mediator, drop = FALSE], kinship = kin, addcovar = covar.med, cores = 4)
      
      
    # Find maximum marker position and ID
    max_marker <- rownames(max_scan1(mediator.qtl, map = map))
     
      
      
    gp[[1]] = gp[[1]][,,max_marker, drop = FALSE]
    # QTL scan on target at mediator's best marker
    target.qtl <- scan1(gp, pheno = overlap.samples[,target.name,drop = FALSE], kinship = kin, addcovar = covar.targ, cores = 4)
      
        
        
        
    targ.pheno = overlap.samples[,target.name,drop = FALSE]
    med.pheno = overlap.samples[,current.mediator,drop = FALSE]
    driver = gp[[1]][,,1]
      
      
    best.MST <- mediation_test(target = targ.pheno,
                               mediator = med.pheno,
                               driver = driver,
                               covar_tar = covar.targ,
                               covar_med = covar.med,
                               kinship = kin,
                               index_name = 'start')
    test <- as.data.frame(best.MST$test)
      
    
    best.marker[i] <- max_marker                                          # Mediator's best cis marker
    sample.size[i] <- nrow(overlap.samples)                               # Sample size between the two
    target.lod[i] <- target.qtl[1,1]                                      # Target LOD score at mediator's best cis marker
    mediator.lod[i] <- mediator.qtl[max_marker,]                          # Mediator LOD score at mediator's best cis marker
    pearson[i] <- cor(overlap.samples)[1,2]                               # Pearson correlation
    spearman[i] <- cor(overlap.samples, method = 'spearman')[1,2]         # Spearman correlation
    best.model[i] <- as.character(test[which.min(test$pvalue),'model'])   # Best model
    best.model.p[i] <-  test[which.min(test$pvalue),'pvalue']             # Best model p-value
    best.triad[i] <- as.character(best.MST$best$triad)                    # Best triad          
    best.triad.p[i] <- best.MST$best$pvalue                               # Best triad p-value
    mediation.lod[i] <- mediation_scan(target = targ.pheno,               # Mediation LOD
                                       mediator = med.pheno,               
                                       driver = driver, 
                                       covar = covar.targ, 
                                       kinship = kin, 
                                       annotation = mediator.annot,
                                       method = 'double-lod-diff', 
                                       index_name = 'start')$lod             
    
    inverse.mediation.lod[i] <- mediation_scan(target = med.pheno,        # Inverse Mediation LOD
                                               mediator = targ.pheno,
                                               driver = driver, 
                                               covar = covar.med,
                                               kinship = kin, 
                                               annotation = target.annot,
                                               method = 'double-lod-diff', 
                                               index_name = 'start')$lod 
    
  
     print(i)
  
}


sparse_data <- cbind(target_and_mediator, best.marker, sample.size, pearson, spearman, target.lod, mediator.lod,
                     mediation.lod, inverse.mediation.lod, best.model,best.model.p, best.triad, best.triad.p)


### Save data as .rds file
saveRDS(sparse_data, file = paste0('attie_islet_protein_mrna_sparse_matrix_chunk_',args[1],'_',args[3],'.rds'))
