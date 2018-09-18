# Define function to pull LOD from linear 


regressionget_lod <- function(y, x, covar){
       n <- nrow(y)
       lm_fit <- lm.fit(y = y, x = as.matrix(cbind(x,covar)))
       rss <- sum(resid(lm_fit)^2)
       tss <- sum((y - mean(y))^2)
       LOD <- (n/2) * log10(tss/rss)
       
       return(LOD)
}

### Load in mRNA and protein data that are in QTL viewer format
load("~/Desktop/Islet Mediation/attie_islet_rZ_qtl_viewer_284.RData")


### Extract required data
mrna <- 'dataset.islet.mrna'
protein <- 'dataset.islet.proteins'

annot.id.protein <- get(protein)$annots %>% dplyr::rename(pos = start)
rownames(annot.id.protein) <- annot.id.protein$protein_id

annot.id.mrna <- get(mrna)$annots %>% dplyr::rename(pos = start)
rownames(annot.id.mrna) <- annot.id.mrna$gene_id

expr.protein <- get(protein)$expr
expr.mrna <- get(mrna)$expr

covar.protein <- get(protein)$covar
covar.mrna <- get(mrna)$covar



### LOD for protein ~ mRNA
protein_id <- annot.id.protein$protein_id
proteins_gene_id <- annot.id.protein$gene_id

protein_on_mrna_lod <- as.data.frame(matrix(NA, nrow = nrow(annot.id.protein), ncol = 1))
colnames(protein_on_mrna_lod) <- c('protein_on_mRNA')

for(i in 1:nrow(annot.id.protein)){
    if(proteins_gene_id[i] %in% annot.id.mrna$gene_id){
       protein_on_mrna_lod[i,] <- regressionget_lod(y = expr.protein[,protein_id[i], drop = FALSE], 
                                                    x = expr.mrna[,proteins_gene_id[i],drop = FALSE], 
                                                    covar = covar.protein)
    }
}
rownames(protein_on_mrna_lod) <- protein_id
saveRDS(protein_on_mrna_lod, file = 'attie_islet_cis_protein_on_mrna_lod_284.rds')




### LOD for mRNA ~ protein
gene_id <- annot.id.mrna$gene_id

mrna_on_protein_lod <- as.data.frame(matrix(NA, nrow = nrow(annot.id.mrna), ncol = 2))
colnames(mrna_on_protein_lod) <- c('protein','mrna_on_protein_lod')

for(i in 1:nrow(annot.id.mrna)){
    if(gene_id[i] %in% annot.id.protein$gene_id){
       med <- annot.id.protein[gene_id[i] == annot.id.protein$gene_id, 'protein_id']
       
       
       temp <- numeric(length = length(med))
       for(j in 1:length(med)){
           temp[j] <- regressionget_lod(y = expr.mrna[,gene_id[i], drop = FALSE], 
                                        x = expr.protein[,med[j],drop = FALSE], 
                                        covar = covar.mrna)
       }
       
      mrna_on_protein_lod[i,] <- cbind(paste0(med, collapse = ','), paste0(temp, collapse = ','))
    }
}
rownames(mrna_on_protein_lod) <- gene_id
saveRDS(mrna_on_protein_lod, file = 'attie_islet_cis_mrna_on_protein_lod_284.rds')


