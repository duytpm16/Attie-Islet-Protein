options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(qtl2)
library(intermediate2)





### Command line arguments/variables to change
args <- commandArgs(trailingOnly = TRUE)
load(args[1])
mediation_data <- readRDS(args[2])
threshold      <- as.numeric(args[3])
target_id      <- args[4]
target_data    <- args[5]
mediator_id    <- args[6]
mediator_data  <- args[7]
chunk_number   <- as.numeric(args[8])
chunk_size     <- as.numeric(args[9])









# Get range of index if chunk_numer and chunk_size is specified
rng <- 1:nrow(mediation_data)

if(!is.null(chunk_number)){
   max_col = nrow(mediation_data)
   rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
  
   if(rng[length(rng)] > max_col) {
      rng = rng[1]:max_col
   }
}


mediation_data <- mediation_data[rng,]











# Get QTL viewer formatted data
annot.id.targ <- get(target_data)$annots
annot.id.targ <- annot.id.targ %>% dplyr::rename(id = target_id) %>% as.data.frame()
annot.id.med  <- get(mediator_data)$annots
annot.id.med  <- annot.id.med %>% dplyr::rename(id = mediator_id) %>% as.data.frame()
covar.targ    <- get(target_data)$covar
covar.med     <- get(mediator_data)$covar
expr.targ     <- get(target_data)$rankz
expr.med      <- get(mediator_data)$rankz










# Find how many rows
n <- nrow(mediation_data)










new_mediation_df <- data.frame()
for(i in 1:n){
  
  
  
  
    # Extract 1 row from mediation_data and turn it into a dataframe by splitting mediator/mediation columns by ','. Subset to those with z-score less than threshold
    temp <- mediation_data[i,]
    mediation_df <- temp %>% 
                          separate_rows(mediator.id, mediator.chr, mediator.start, mediator.end, mediator.symbol, mediation.lod, sep = ',') %>%
                          mutate(mediation.lod = as.numeric(mediation.lod),
                                 mediator.start = as.numeric(mediator.start),
                                 mediator.end = as.numeric(mediator.end),
                                 mediation.z.lod = scale(mediation.lod)) %>%
                          subset(mediator.id == annot.id.targ[temp$target.id, mediator_id])
  
  
  
  
  
  
  
    # Do inverse mediation and qtl2 scan on those that pass the threshold and finally add to new_mediation_df
    if(nrow(mediation_df) != 0){
       gp = genoprobs[,mediation_df$target.qtl.chr[1]]
       gp[[1]] = gp[[1]][,,mediation_df$marker.id[1], drop = FALSE]
       temp.expr.targ <- expr.targ[, mediation_df$target.id[1], drop = FALSE]
       for(j in 1:nrow(mediation_df)){
        
           # Inverse mediation
           mediation_df$inv.mediation.lod[j] <- mediation_scan(target     = expr.med[, mediation_df$mediator.id[j],drop = FALSE],
                                                               mediator   = temp.expr.targ,
                                                               annotation = annot.id.targ[mediation_df$target.id[j],],
                                                               driver     = genoprobs[[mediation_df$target.qtl.chr[j]]][,, mediation_df$marker.id[j]],
                                                               kinship    = K[[mediation_df$target.qtl.chr[j]]],
                                                               covar      = covar.med,
                                                               method     = 'double-lod-diff',
                                                               verbose    = FALSE,
                                                               index_name = 'start')$lod
        
        
           # qtl2 scan1
           intersect.samples <- intersect(rownames(temp.expr.targ[complete.cases(temp.expr.targ),,drop = FALSE]), rownames(expr.med[, mediation_df$mediator.id[j],drop = FALSE]))
           mediation_df$mediator.lod[j] <- scan1(pheno = expr.med[,mediation_df$mediator.id[j],drop = FALSE],
                                                 genoprobs = gp,
                                                 kinship = K[[mediation_df$target.qtl.chr[j]]],
                                                 addcovar = covar.med)
      }     
      
    new_mediation_df <- bind_rows(new_mediation_df, mediation_df)
    }
    
    
    print(paste0(i,'of',n))
}




# Calculate target and mediator lod drop and lod drop proportion
new_mediation_df$target.drop   <- new_mediation_df$target.lod - new_mediation_df$mediation.lod
new_mediation_df$target.prop   <- new_mediation_df$target.drop / new_mediation_df$target.lod
new_mediation_df$mediator.drop <- new_mediation_df$mediator.lod - new_mediation_df$inv.mediation.lod
new_mediation_df$mediator.prop <- new_mediation_df$mediator.drop / new_mediation_df$mediator.lod





### Rearrange column of data
new_mediation_df <- new_mediation_df %>% select(target.id, target.chr, target.start, target.end, target.symbol, 
                                                target.qtl.chr,target.qtl.pos, marker.id, target.cis,
                                                mediator.id, mediator.chr, mediator.start, mediator.end, mediator.symbol, mediation.z.lod,
                                                target.lod, mediator.lod, mediation.lod, inv.mediation.lod,
                                                target.drop, target.prop, mediator.drop, mediator.prop)



### Extract only cis QTLs
new_mediation_df <- new_mediation_df[new_mediation_df$target.cis, ]





### Save data as .rds
file_name <- c(target_id, mediator_id)      # Just for naming conventions
file_name <- gsub('_id','',file_name)
file_name <- gsub('gene','mrna', file_name)


if(!is.null(chunk_number)){
   saveRDS(new_mediation_df, paste0('attie_islet_',file_name[1],'_filter_cis_yandell_intermediate_',file_name[2],'_chunk_',chunk_number,'.rds'))
}else{
   saveRDS(new_mediation_df, paste0('attie_islet_',file_name[1],'_filter_cis_yandell_intermediate_',file_name[2],'.rds'))
}
