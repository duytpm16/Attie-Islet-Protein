options(stringsAsFactors = FALSE)



library(crayon, lib.loc = '~/Rlibs')
library(dplyr, lib.loc = '~/Rlibs')
library(data.table, lib.loc = '~/Rlibs')
library(intermediate, lib.loc = '~/Rlibs')

load('attie_islet_284_qtl_viewer.RData')


### Extract required data from QTL viewer dataset
intermediate_function <- function(target_id, target_data, mediator_data, mediator_id = NULL, chunk_number = NULL, chunk_size = NULL){

       if(target_data == mediator_data){
          
          # Get data
          annot.id <- get(target_data)$annots
          annot.id <- annot.id %>% dplyr::rename(pos = start) %>% as.data.frame()
	  rownames(annot.id) <- annot.id[,target_id]
          covar <- get(target_data)$covar
          expr <- get(target_data)$rankz
          
          
          lod.peaks <- get(target_data)$lod.peaks$additive
          
       
          
          # Get range of index
          rng <- 1:nrow(lod.peaks)
          
          if(!is.null(chunk_number)){
            max_col = nrow(lod.peaks)
            rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
            if(rng[length(rng)] > max_col) {
              
              rng = rng[1]:max_col
              
            }
          }
         


          # Get the distal pQTL lodpeaks info
          qtl.annot <- lod.peaks$annot.id
          qtl.chr <- lod.peaks$qtl.chr
          qtl.marker <- lod.peaks$marker.id
          annot.index <- match(qtl.annot, annot.id[,target_id])
          expr.index <- match(qtl.annot, colnames(expr))
          
          
          mediation_results <- list()
          for(i in rng){
        
        
              # Mediation: Y ~ Q  + M + covar
              med <- mediation.scan(target = expr[, qtl.annot[i], drop = FALSE],
                                    mediator = expr[,-expr.index[i]],
                                    annotation = annot.id[-annot.index[i],],
                                    qtl.geno = genoprobs[[qtl.chr[i]]][,,qtl.marker[i]],
                                    covar = covar,
                                    method = "double-lod-diff",
                                    verbose = FALSE)
      
      
              # Save results
              result <- data.frame(target = qtl.annot[i], 
                                   target.qtl.chr = qtl.chr[i],
                                   target.qtl.pos = lod.peaks$qtl.pos[i],
			           marker.id = lod.peaks$marker.id[i],
                                   target.lod = lod.peaks$lod[i],
                                   target.chr = lod.peaks$gene.chr[i],
                                   target.start = lod.peaks$gene.start[i],
                                   target.end = lod.peaks$gene.end[i],
                                   target.symbol = lod.peaks$gene.symbol[i],
                                   target.cis = lod.peaks$cis[i],
                                   mediator = paste0(med[,mediator_id], collapse = ','),
                                   mediator.chr = paste0(med$chr, collapse = ','),
                                   mediator.start = paste0(med$pos, collapse = ','),
                                   mediator.end = paste0(med$end, collapse = ','),
                                   mediator.symbol = paste0(med$symbol, collapse = ','),
                                   mediation.lod = paste0(med$LOD, collapse = ','))
                                   
              
              mediation_results[[i]] <- result
              print(paste(i,'out of', nrow(lod.peaks)))
         }
         
          mediation_results <- rbindlist(mediation_results)
          
          
          ### Save data as .rds
          if(target_id == 'protein_id'){
            if(!is.null(chunk_number)){
              saveRDS(mediation_results, paste0('attie_islet_protein_intermediate_protein_chunk_',chunk_number,'.rds'))
            }else{
              saveRDS(mediation_results, 'attie_islet_protein_intermediate_protein.rds')
            }
          }
          
          if(target_id == 'gene_id'){
            if(!is.null(chunk_number)){
              saveRDS(mediation_results, paste0('attie_islet_mrna_intermediate_mrna_chunk_',chunk_number,'.rds'))
            }else{  
              saveRDS(mediation_results, 'attie_islet_mrna_intermediate_mrna.rds')
            }
          }
      } # if target_data == mediator_data
  
  
  
  
  
  
  
  
  
  
  
    
    if(target_data != mediator_data){
       
       # Get data
       annot.id.targ <- get(target_data)$annots
       annot.id.targ <- annot.id.targ %>% dplyr::rename(pos = start) %>% as.data.frame()
       annot.id.med <- get(mediator_data)$annots
       annot.id.med <- annot.id.med %>% dplyr::rename(pos = start) %>% as.data.frame()
       rownames(annot.id.med) <- annot.id.med[,mediator_id]
       covar.targ <- get(target_data)$covar
       covar.med <- get(mediator_data)$covar
       expr.targ <- get(target_data)$rankz
       expr.med <- get(mediator_data)$rankz
       lod.peaks <- get(target_data)$lod.peaks$additive
       
       

       
       # Get range of index
       rng <- 1:nrow(lod.peaks)
       
       if(!is.null(chunk_number)){
         max_col = nrow(lod.peaks)
         rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
         if(rng[length(rng)] > max_col) {
           
           rng = rng[1]:max_col
           
         }
       }
       
       
       qtl.annot <- lod.peaks$annot.id
       qtl.chr <- lod.peaks$qtl.chr
       qtl.marker <- lod.peaks$marker.id
       mediation_results <- list()
       for(i in rng){
           # Mediation: Y ~ Q  + M + covar
           med <- mediation.scan(target = expr.targ[, qtl.annot[i], drop = FALSE],
                                 mediator = expr.med,
                                 annotation = annot.id.med,
                                 qtl.geno = genoprobs[[qtl.chr[i]]][,,qtl.marker[i]],
                                 covar = covar.targ,
        		         method = "double-lod-diff",
                                 verbose = FALSE)
           
           # Save results
           result <- data.frame(target = qtl.annot[i], 
                                target.qtl.chr = qtl.chr[i],
                                target.qtl.pos = lod.peaks$qtl.pos[i],
                                marker.id = lod.peaks$marker.id[i],
				target.lod = lod.peaks$lod[i],
                                target.chr = lod.peaks$gene.chr[i],
                                target.start = lod.peaks$gene.start[i],
                                target.end = lod.peaks$gene.end[i],
                                target.symbol = lod.peaks$gene.symbol[i],
                                target.cis = lod.peaks$cis[i],
                                mediator = paste0(med[,mediator_id], collapse = ','),
                                mediator.chr = paste0(med$chr, collapse = ','),
                                mediator.start = paste0(med$pos, collapse = ','),
                                mediator.end = paste0(med$end, collapse = ','),
                                mediator.symbol = paste0(med$symbol, collapse = ','),
                                mediation.lod = paste0(med$LOD, collapse = ','))
         
           mediation_results[[i]] <- result
           print(paste0(i,'out of', nrow(lod.peaks)))
       }
       
       
       mediation_results <- rbindlist(mediation_results)
       
       ### Save data as .rds
       if(target_id == 'protein_id'){
         if(!is.null(chunk_number)){
           saveRDS(mediation_results, paste0('attie_islet_protein_intermediate_mrna_chunk_',chunk_number,'.rds'))
         }else{
           saveRDS(mediation_results, 'attie_islet_protein_intermediate_mrna.rds')
         }
       }
       
       if(target_id == 'gene_id'){
         if(!is.null(chunk_number)){
           saveRDS(mediation_results, paste0('attie_islet_mrna_intermediate_protein_chunk_',chunk_number,'.rds'))
         }else{  
           saveRDS(mediation_results, 'attie_islet_mrna_intermediate_protein.rds')
         }
       }
   } # if target_data != mediator_data
}      


args = commandArgs(trailingOnly = TRUE)
intermediate_function(args[1],args[2],args[3],args[4], as.numeric(args[5]), as.numeric(args[6]))



