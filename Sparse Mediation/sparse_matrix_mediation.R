options(stringsAsFactors = FALSE)
library(dplyr)
library(intermediate2)



### Arguments to pass
args = commandArgs(trailingOnly = TRUE)
chunk_number = as.numeric(args[1])


### Read in the sparse matrix
sparse_matrix <- readRDS(paste0('attie_islet_sparse_matrix_',args[1],'.rds'))



### Extract rankZ, annotation, and covariate data
rankz.protein <- dataset.islet.proteins$rankz
annot.protein <- dataset.islet.proteins$annots
covar.protein <- dataset.islet.proteins$covar




### Get target and mediator ids
annot.protein <- annot.protein %>% dplyr::rename(id = protein_id)

target.id <- sparse_matrix$target.id

mediator.id <- sparse_matrix$mediator.id
mediator.chr <- sparse_matrix$chr
mediator.marker <- sparse_matrix$marker.id




### Create new vectors to store mediation analysis results
n <- nrow(sparse_matrix)
best.model <- character(length = n)
best.p <- numeric(length = n)
best.triad <- character(length = n)
best.triad.p <- numeric(length = n)
mediation.lod <- numeric(length = n)
inverse.mediation.lod <- numeric(length = n)






### Mediation Begin
for(i in 1:nrow(sparse_matrix)){
        
        # Extracting data
        targ.pheno = rankz.protein[,target.id[i], drop = FALSE]
        med.pheno = rankz.protein[,mediator.id[i], drop = FALSE]
        kin = K[[mediator.chr[i]]]
        targ.annot <- annot.protein[target.id[i],]
        med.annot <- annot.protein[mediator.id[i],]
        driver = genoprobs[[mediator.chr[i]]][,,mediator.marker[i]]
        
        
        
        
        
        # Find best causalMST method: best overall model and pvalue, and best triad and triad pvalue (no undecided)
        best.MST <- mediation_test(target = targ.pheno,
                                   mediator = med.pheno,
                                   driver = driver,
                                   covar_tar = covar.protein,
                                   covar_med = covar.protein,
                                   kinship = kin,
                                   index_name = 'start')
        test <- as.data.frame(best.MST$test)
        
        best.model[i] <- as.character(test[which.min(test$pvalue),'model'])
        best.p[i] <- test[which.min(test$pvalue),'pvalue']
        best.triad[i] <- as.character(best.MST$best$triad)
        best.triad.p[i] <- best.MST$best$pvalue
        
        
        
        # Mediation
        mediation.lod[i] <- mediation_scan(target = targ.pheno,
                                           mediator = med.pheno,
                                           driver = driver,
                                           covar = covar.protein,
                                           kinship = kin,
                                           annotation = med.annot,
                                           method = 'double-lod-diff',
                                           index_name = 'start')$lod
        
        # Inverse Mediation
        inverse.mediation.lod[i] <- mediation_scan(target = med.pheno,
                                                   mediator = targ.pheno,
                                                   driver = driver,
                                                   covar = covar.protein,
                                                   kinship = kin,
                                                   annotation = targ.annot,
                                                   method = 'double-lod-diff',
                                                   index_name = 'start')$lod
        
        
        print(i)
    
}

sparse_matrix <- cbind(sparse_matrix, mediation.lod, inverse.mediation.lod, best.model, best.p, best.triad, best.triad.p)
saveRDS(sparse_matrix, file = paste0('attie_islet_protein_sparse_matrix_chunk_',args[1],'.rds'))

