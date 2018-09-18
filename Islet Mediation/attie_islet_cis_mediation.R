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



### Protein cis mediation and inverse mediation with mRNA
# 
#   Extracting some columns from protein annots so the code doesn't look messy
protein_id <- annot.id.protein$protein_id
proteins_gene_id <- annot.id.protein$gene_id
chr = annot.id.protein$chr
nearest.marker <- annot.id.protein$nearest.marker.id
ids <- annot.id.protein[,'protein_id',drop = TRUE]

protein_cis_med_lod <- as.data.frame(matrix(NA, nrow = nrow(annot.id.protein), ncol = 2))
colnames(protein_cis_med_lod) <- c('cis_mediator_mRNA','cis_inv_mediation_mRNA')

for(i in 1:nrow(annot.id.protein)){
  
    if(proteins_gene_id[i] %in% annot.id.mrna$gene_id){
       med <- proteins_gene_id[i]
       temp <- mediation.scan(target = expr.protein[, protein_id[i], drop = FALSE],
                              mediator = expr.mrna[, med, drop = FALSE],
                              annotation = annot.id.mrna[proteins_gene_id[i],],
                              qtl.geno = genoprobs[[chr[i]]][,,nearest.marker[i]],
                              covar = covar.protein)
       
       
       temp.inv <- mediation.scan(target = expr.mrna[, proteins_gene_id[i], drop = FALSE],
                                  mediator = expr.protein[, protein_id[i], drop = FALSE],
                                  annotation = annot.id.protein[protein_id[i],],
                                  qtl.geno = genoprobs[[chr[i]]][,,nearest.marker[i]],
                                  covar = covar.mrna)
       
       protein_cis_med_lod[i,] <- cbind(temp[,'LOD'],temp.inv[,'LOD'])

    }
}
rownames(protein_cis_med_lod) <- protein_id
saveRDS(protein_cis_med_lod, file = 'attie_islet_protein_cis_mediation_lod_284.rds')





### mRNA cis mediation and inverse mediation with mrna
# 
#   Extracting some columns from mrna annots so the code doesn't look messy
gene_id <- annot.id.mrna$gene_id
chr = annot.id.mrna$chr
nearest.marker <- annot.id.mrna$nearest.marker.id
ids <- annot.id.mrna[,'gene_id',drop = TRUE]

mrna_cis_med_lod <- as.data.frame(matrix(NA, nrow = nrow(annot.id.mrna), ncol = 3))
colnames(mrna_cis_med_lod) <- c('mediator_protein','cis_mediator_protein','cis_inv_mediation_protein')

for(i in 1:nrow(annot.id.mrna)){
  
    if(gene_id[i] %in% annot.id.protein$gene_id){
        med <- annot.id.protein[which(gene_id[i] == annot.id.protein$gene_id), 'protein_id']
        temp <- mediation.scan(target = expr.mrna[, gene_id[i], drop = FALSE],
                               mediator = expr.protein[, med, drop = FALSE],
                               annotation = annot.id.protein[med,],
                               qtl.geno = genoprobs[[chr[i]]][,,nearest.marker[i]],
                               covar = covar.mrna)
      
      
        temp.inv <- as.data.frame(matrix(NA, nrow = length(med), ncol = 10))
        for(j in 1:length(med)){
            temp.inv[j,] <- mediation.scan(target = expr.protein[, med[j], drop = FALSE],
                                       mediator = expr.mrna[, gene_id[i], drop = FALSE],
                                       annotation = annot.id.mrna[gene_id[i],],
                                       qtl.geno = genoprobs[[chr[i]]][,,nearest.marker[i]],
                                       covar = covar.protein)
        }
      
        mrna_cis_med_lod[i,] <- cbind(paste0(med, collapse = ','), paste0(temp[,'LOD'],collapse = ','), paste0(temp.inv[,ncol(temp.inv)],collapse = ','))
    }
  
  
}
rownames(mrna_cis_med_lod) <- gene_id
saveRDS(mrna_cis_med_lod, file = 'attie_islet_mrna_cis_mediation_lod_284.rds')

    
