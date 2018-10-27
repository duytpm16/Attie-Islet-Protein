options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(intermediate)
library(qtl2)




### Command line arguments/variables to change
# args <- commandArgs(trailingOnly = TRUE)
# load(args[1])
# mediation_data <- args[2]
# threshold      <- args[3]
# target_id      <- args[4]
# target_data    <- args[5]
# mediator_id    <- args[6]
# mediator_data  <- args[7]








  
# Get QTL viewer formatted data
annot.id.targ <- get(target_data)$annots
annot.id.targ <- annot.id.targ %>% dplyr::rename(pos = start) %>% as.data.frame()
annot.id.med  <- get(mediator_data)$annots
annot.id.med  <- annot.id.med %>% dplyr::rename(pos = start) %>% as.data.frame()
covar.targ    <- get(target_data)$covar
covar.med     <- get(mediator_data)$covar
expr.targ     <- get(target_data)$rankz
expr.med      <- get(mediator_data)$rankz

  
  

stopifnot(colnames(expr.targ) == annot.id.targ[, target_id])
stopifnot(colnames(expr.med) == annot.id.med[, mediator_id])




# Find how many rows
n <- nrow(mediation_data)
  
  










# Extract 1 row from mediation_data and turn it into a dataframe by splitting some columns by ,
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
                        subset(target.qtl.chr == mediator.chr & abs(target.qtl.pos - mediator.start) <= 10 & (abs(mediation.z.lod) > abs(threshold)))
    
    
    
    
    
    
    
    # Do inverse mediation and qtl2 scan on those that pass the threshold and finally add to new_mediation_df
    if(nrow(mediation_df) != 0){

       for(j in 1:nrow(mediation_df)){
           
           # Inverse mediation
           mediation_df$inv.mediation.lod[j] <- mediation.scan(target = expr.med[,mediation_df$mediator.id[j],drop = FALSE],
                                                               mediator = expr.targ[,mediation_df$target.id[j], drop = FALSE],
                                                               annotation = annot.id.targ[mediation_df$target.id[j],],
                                                               qtl.geno = genoprobs[[mediation_df$target.qtl.chr[j]]][,,mediation_df$marker.id[j]],
                                                               covar = covar.med,
                                                               method = 'double-lod-diff')$LOD
           
           
           # qtl2 scan1
           gp = genoprobs[,mediation_df$target.qtl.chr[j]]
           gp[[1]] = gp[[1]][,,mediation_df$marker.id[j], drop = FALSE]
           mediation_df$mediator.lod[j] <- scan1(pheno = expr.med[,mediation_df$mediator.id[j],drop = FALSE],
                                                 genoprobs = gp,
                                                 kinship = K[[mediation_df$target.qtl.chr[j]]],
                                                 addcovar = covar.med)
       }     
      
       
       # Calculate target and mediator lod drop and lod drop proportion
       mediation_df$target.drop   <- mediation_df$target.lod - mediation_df$mediation.lod
       mediation_df$target.prop   <- mediation_df$target.drop / mediation_df$target.lod
       mediation_df$mediator.drop <- mediation_df$mediator.lod - mediation_df$inv.mediation.lod
       mediation_df$mediator.prop <- mediation_df$mediator.drop / mediation_df$mediator.lod
      
       new_mediation_df <- bind_rows(new_mediation_df, mediation_df)
   }
    
    
   print(paste0(i,'of',n))
}









### Rearrange column of data
new_mediation_df <- new_mediation_df %>% select(target.id, target.chr, target.start, target.end, target.symbol, 
                                                target.qtl.chr,target.qtl.pos, marker.id, target.cis,
                                                mediator.id, mediator.chr, mediator.start, mediator.end, mediator.symbol, mediation.z.lod,
                                                target.lod, mediator.lod, mediation.lod, inv.mediation.lod,
                                                target.drop, target.prop, mediator.drop, mediator.prop)









### Save data as .rds
file_name <- c(target_id, mediator_id)      # Just for naming conventions
file_name <- gsub('_id','',file_name)
file_name <- gsub('gene','mrna', file_name)


saveRDS(new_mediation_df, paste0('attie_islet_',file_name[1],'_filtered_simecek_intermediate_',file_name[2],'.rds'))

