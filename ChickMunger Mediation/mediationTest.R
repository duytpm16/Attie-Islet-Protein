options(stringsAsFactors = FALSE)
options(scipen = 999)

library(crayon, lib.loc = '~/Rlibs')
library(bindrcpp, lib.loc = '~/Rlibs')
library(intermediate2, lib.loc = '~/Rlibs')
library(dplyr, lib.loc = '~/Rlibs')
library(data.table, lib.loc = '~/Rlibs')


load('attie_islet_284_qtl_viewer.RData')


### Extract required data from QTL viewer dataset
mediationTest_function <- function(target_id, target_data, mediator_data, mediator_id, mediation_table, threshold = NULL, chunk_number = NULL, chunk_size = NULL){
  
    
     # Get data
     annot.id.targ <- get(target_data)$annots
     annot.id.targ <- annot.id.targ %>% dplyr::rename(pos = start) %>% as.data.frame()
     colnames(annot.id.targ)[colnames(annot.id.targ) == target_id] <- 'id'
     annot.id.med <- get(mediator_data)$annots
     annot.id.med <- annot.id.med %>% dplyr::rename(pos = start) %>% as.data.frame()
     colnames(annot.id.med)[colnames(annot.id.med) == mediator_id] <- 'id'
     covar.targ <- get(target_data)$covar
     covar.med <- get(mediator_data)$covar
     expr.targ <- get(target_data)$rankz
     expr.med <- get(mediator_data)$rankz
    
    
    
    
     # Get range of index
     rng <- 1:nrow(mediation_table)
    
     if(!is.null(chunk_number)){
        max_col = nrow(mediation_table)
        rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
        if(rng[length(rng)] > max_col) {
        
           rng = rng[1]:max_col
        
        }
      }
      mediation_table <- mediation_table[rng,] 
    
      target <- mediation_table$target
      mediator <- mediation_table$mediator
      target.qtl.chr <- mediation_table$target.qtl.chr
      target.qtl.marker <- mediation_table$marker.id
    
      mediation_results <- list()
      for(i in 1:nrow(mediation_table)){
      

      
          ### Mediation: Y ~ Q  + M + covar
          med <- mediation_test(target = expr.targ[, target[i], drop = FALSE],
                                mediator = expr.med[,mediator[i], drop = FALSE],
                                annotation = annot.id.med[annot.id.med$id %in% mediator[i], ],
                                driver = genoprobs[[target.qtl.chr[i]]][,,target.qtl.marker[i]],
                                covar_tar = covar.targ,
                                covar_med = covar.med,
                                verbose = FALSE,
				minN = 0)
      
          med.best <- as.data.frame(med$test)
          ### Save results
          result <- data.frame(p.causal = med$test[med$test$model %in% 'causal', 'pvalue', drop = TRUE],
                               p.reactive = med$test[med$test$model %in% 'reactive', 'pvalue', drop = TRUE],
                               p.independent = med$test[med$test$model %in% 'independent', 'pvalue', drop = TRUE],
                               p.undecided = med$test[med$test$model %in% 'undecided', 'pvalue', drop = TRUE],
                               best.model = as.character(med.best[which.min(med.best$pvalue),'model']),
                               p.best.model = as.character(med.best[which.min(med.best$pvalue),'pvalue']),
                               best.mediation = med$best$mediation,
                               best.LRmed = med$best$LRmed,
                               best.LR = med$best$LR,
                               t.d_t = med$fitsLR$t.d_t,
                               m.d_m = med$fitsLR$m.d_m,
                               t.m_t = med$fitsLR$t.m_t,
                               m.t_m = med$fitsLR$m.t_m,
                               t.md_t.m = med$fitsLR$t.md_t.m,
                               t.md_t = med$fitsLR$t.md_t)
        
          mediation_results[[i]] <- result
          print(paste0(i,'out of', nrow(mediation_table)))
     }
    
    
     mediation_results <- rbindlist(mediation_results)
     mediation_results <- cbind(mediation_table, mediation_results)
     return(mediation_results)
  
}

args = commandArgs(trailingOnly = TRUE)

mediation_table = readRDS(args[5])



### Protein-Protein MST
#result_protein_protein <- mediationTest_function(target_id = args[1], target_data = args[2], 
#                                                 mediator_data = args[3], mediator_id = args[4],
#                                                 mediation_table = mediation_table)
#   Save protein-protein MST
#saveRDS(result_protein_protein, file = 'attie_islet_protein_protein_284_casualMST.rds')
#write.table(result_protein_protein, file = 'attie_islet_protein_protein_284_casualMST.txt', sep = '\t', row.names = FALSE, col.names = TRUE)




### Protein-mRNA MST
result <- mediationTest_function(target_id = args[1], target_data = args[2],
                                              mediator_data = args[3], mediator_id = args[4],
                                              mediation_table = mediation_table, chunk_number = as.numeric(args[6]), chunk_size = as.numeric(args[7]))


#   Save protein-mrna MST
saveRDS(result, file = paste0('attie_islet_protein_mrna_284_abs_4_causalMST_chunk_',args[6],'.rds'))
#write.table(result, file = paste0('attie_islet_',targ,'_',med,'_284_abs_4_causalMST.txt'), sep = '\t', row.names = FALSE, col.names = TRUE)




### mRNA-Protein MST
#result_mrna_protein <- mediationTest_function(target_id = args[1], target_data = args[2],
#                                              mediator_data = args[3], mediator_id = args[4],
#                                              mediation_table = mediation_table)
#   Save mRNA-Protein MST
#saveRDS(result_mrna_protein, file = 'attie_islet_mrna_protein_284_casualMST.rds')
#write.table(result_mrna_protein, file = 'attie_islet_mrna_protein_284_casualMST.txt', sep = '\t', row.names = FALSE, col.names = TRUE)




### mRNA-mRNA MST
#result_mrna_mrna <- mediationTest_function(target_id = args[1], target_data = args[2],
#                                           mediator_data = args[3], mediator_id = args[4],
#                                           mediation_table = mediation_table)
#   Save mRNA-mRNA MST
#saveRDS(result_mrna_mrna, file = 'attie_islet_mrna_mrna_284_casualMST.rds')
#write.table(result_mrna_mrna, file = 'attie_islet_mrna_mrna_284_casualMST.txt', sep = '\t', row.names = FALSE, col.names = TRUE)




