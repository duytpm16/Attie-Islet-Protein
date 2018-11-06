######################################################################################################################
#
#   This script filter's the mediation dataframe generated from simecek_mediation_scan.R given a threshold. Qtl2 scans, 
#        inverse mediation, LOD drop, and LOD drop proportion is performed on all mediator with global 
#        mediation Z lod score < threshold.
#
#   Petr Simecek's intermediate package: install_github('simecek/intermediate')
#
#
#   Input:
#      1.) Load in the QTL Viewer environment
#      2.) mediation_data : dataframe generated from yandell_mediation_scan.R
#      3.) threshold      : cut off theshold to filter mediators
#      4.) target_id      : Name of the column in annots that will be used to select the targets
#      5.) target_data    : dataset.* list where annots contain the column name in [2]
#      6.) mediator_id    : Name of the column in annots that will be used to select the mediators
#      7.) mediator_data  : dataset.* list where annots contain the column name in [4]
#      8.) chunk_number   : numeric value indicating which portion of the lod peaks table to break. Need a fixed chunk_size
#      9.) chunk_size     : numeric value indicating the size to break the lod peaks table
#
#
#
#   Output:
#      1.) Long dataframe version of the mediation_data for all mediators that pass threshold. 
#         
#
#
#
#   Author: Duy Pham
#   Date:   November 6, 2018
#   E-mail: duy.pham@jax.org
#
######################################################################################################################
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(intermediate)
library(qtl2)




### Command line arguments/variables to change
args <- commandArgs(trailingOnly = TRUE)
load(args[1])                             # Load QTL Viewer .RData file
mediation_data <- readRDS(args[2])        # Load mediation .rds file
threshold      <- as.numeric(args[3])     # Z-score threshold
target_id      <- args[4]                 # Either protein_id/gene_id
target_data    <- args[5]                 # Which dataset.* to use as target
mediator_id    <- args[6]                 # Either protein_id/gene_id
mediator_data  <- args[7]                 # Which dataset.* to use as mediator
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
                        subset(target.qtl.chr == mediator.chr & abs(target.qtl.pos - mediator.start) <= 10 & mediation.z.lod < threshold)
    
    
    
    
    
    
    
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









### Save data as .rds
### Save data as .rds
file_name <- c(target_id, mediator_id)      # Just for naming conventions
file_name <- gsub('_id','',file_name)
file_name <- gsub('gene','mrna', file_name)


if(!is.null(chunk_number)){
   saveRDS(new_mediation_df, paste0('attie_islet_',file_name[1],'_filter_simecek_intermediate_',file_name[2],'_chunk_',chunk_number,'.rds'))
}else{
   saveRDS(new_mediation_df, paste0('attie_islet_',file_name[1],'_filter_simecek_intermediate_',file_name[2],'.rds'))
}
