options(stringsAsFactors = FALSE)
library(dplyr)
library(intermediate2)



### Variables to change
args = commandArgs(trailingOnly = TRUE)
load("~/Desktop/Attie Final/Islet Mediation/attie_islet_284_qtl_viewer.RData")
sparse_matrix <- readRDS(paste0('attie_islet_sparse_matrix_1000.rds'))
chunk_number = as.numeric(1)
chunk_size = as.numeric(1000)
targ <- 'dataset.islet.proteins'
med <- 'dataset.islet.proteins'
targ_id <- 'protein_id'
med_id <- 'protein_id'




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



# Get chunk range for target
max_col = nrow(sparse_matrix)
rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
if(rng[length(rng)] > max_col) {
   rng = rng[1]:max_col
}




sparse_matrix <- sparse_matrix[rng,]


### Get target and mediator ids
target.id <- sparse_matrix$target.id

mediator.id <- sparse_matrix$mediator.id
mediator.chr <- sparse_matrix$chr
mediator.marker <- sparse_matrix$best.marker




### Create new vectors to store mediation analysis results
n <- nrow(sparse_matrix)
best.model <- character(length = n)
best.p <- numeric(length = n)
best.triad <- character(length = n)
best.triad.p <- numeric(length = n)
causal.p <- numeric(length = n)
reactive.p <- numeric(length = n)
independent.p <- numeric(length = n)
undecided.p <- numeric(length = n)
t.d_t <- numeric(length = n)
m.d_m <- numeric(length = n)
t.m_t <- numeric(length = n)
m.t_m <- numeric(length = n)
t.md_t.m <- numeric(length = n)
t.md_t <- numeric(length = n)
mediation.lod <- numeric(length = n)
inverse.mediation.lod <- numeric(length = n)







for(i in 1:nrow(sparse_matrix)){
        
        targ.pheno = rankz.targ[,target.id[i], drop = FALSE]
        med.pheno = rankz.med[,mediator.id[i], drop = FALSE]
        kin = K[[mediator.chr[i]]]
        targ.annot <- annot.targ[target.id[i],]
        med.annot <- annot.med[mediator.id[i],]
        driver = genoprobs[[mediator.chr[i]]][,,mediator.marker[i]]
        
        
        
        
        ### Mediation Analyis
        # Get mediation_test results
        best.MST <- mediation_test(target = targ.pheno,
                                   mediator = med.pheno,
                                   driver = driver,
                                   covar_tar = covar.targ,
                                   covar_med = covar.med,
                                   kinship = kin,
                                   index_name = 'start')
        test <- as.data.frame(best.MST$test)
        
        
        
        
        
        
        
        
        ### Store info
        best.model[i] <- as.character(test[which.min(test$pvalue),'model'])
        best.p[i] <- test[which.min(test$pvalue),'pvalue']
        best.triad[i] <- as.character(best.MST$best$triad)
        best.triad.p[i] <- best.MST$best$pvalue
        causal.p[i] <- test[test$model == 'causal', 'pvalue']
        reactive.p[i] <- test[test$model == 'reactive', 'pvalue']
        independent.p[i] <- test[test$model == 'independent', 'pvalue']
        undecided.p[i] <- test[test$model == 'undecided', 'pvalue']
        t.d_t[i] <- best.MST$fitsLR$t.d_t
        m.d_m[i] <- best.MST$fitsLR$m.d_m
        t.m_t[i] <- best.MST$fitsLR$t.m_t
        m.t_m[i] <- best.MST$fitsLR$m.t_m
        t.md_t.m[i] <- best.MST$fitsLR$t.md_t.m
        t.md_t[i] <- best.MST$fitsLR$t.md_t
        
        
        
        
        
        # Mediation
        mediation.lod[i] <- mediation_scan(target = targ.pheno,
                                           mediator = med.pheno,
                                           driver = driver,
                                           covar = covar.targ,
                                           kinship = kin,
                                           annotation = med.annot,
                                           method = 'double-lod-diff',
                                           index_name = 'start')$lod
        
        # Inverse Mediation
        inverse.mediation.lod[i] <- mediation_scan(target = med.pheno,
                                                   mediator = targ.pheno,
                                                   driver = driver,
                                                   covar = covar.med,
                                                   kinship = kin,
                                                   annotation = targ.annot,
                                                   method = 'double-lod-diff',
                                                   index_name = 'start')$lod
        
        
        print(i)
    
}

sparse_matrix <- cbind(sparse_matrix, mediation.lod, inverse.mediation.lod, best.model, best.p, best.triad, best.triad.p, causal.p, reactive.p, independent.p,
                       undecided.p, t.d_t, m.d_m, t.m_t, m.t_m, t.md_t.m, t.md_t)
saveRDS(sparse_matrix, file = paste0('attie_islet_protein_sparse_matrix_complete',args[2],'.rds'))

