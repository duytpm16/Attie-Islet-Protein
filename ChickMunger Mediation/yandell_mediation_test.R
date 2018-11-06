######################################################################################################################
#
#   This script performs causalMST with a kinship matrix given two QTL Viewer formatted datasets using
#        Bryan Yandell's intermediate package: install_github('byandell/intermediate')
#
#
#   Input:
#      1.) Load in the QTL Viewer environment
#      2.) mediation_data :  Filtered mediation data generated from one of the filter mediation scripts
#      3.) target_id      :  Name of the column in annots that will be used to select the targets
#      4.) target_data    :  dataset.* list where annots contain the column name in [2]
#      5.) mediator_id    :  Name of the column in annots that will be used to select the mediators
#      6.) mediator_data  :  dataset.* list where annots contain the column name in [4]
#      7.) chunk_number   :  numeric value indicating which portion of the lod peaks table to break. Need a fixed chunk_size
#      8.) chunk_size     :  numeric value indicating the size to break the lod peaks table
#
#
#
#   Output:
#      1.) Dataframe with causalMST results using the mediation_test function
#
#
#
#   Author: Duy Pham
#   Date:   November 6, 2018
#   E-mail: duy.pham@jax.org
#
######################################################################################################################

### Load library packages
options(stringsAsFactors = FALSE)
options(scipen = 999)


library(intermediate2)
library(dplyr)







### Command line arguments/variables to change
args <- commandArgs(trailingOnly = TRUE)
load(args[1])
mediation_data <- readRDS(args[2])
target_id      <- args[3]
target_data    <- args[4]
mediator_id    <- args[5]
mediator_data  <- args[6]
chunk_number   <- as.numeric(args[7])
chunk_size     <- as.numeric(args[8])






# Get data for mediation_test function
annot.id.targ <- get(target_data)$annots
annot.id.targ <- annot.id.targ %>% dplyr::rename(id = target_id) %>% as.data.frame()
annot.id.med  <- get(mediator_data)$annots
annot.id.med  <- annot.id.med %>% dplyr::rename(id = mediator_id) %>% as.data.frame()
covar.targ    <- get(target_data)$covar
covar.med     <- get(mediator_data)$covar
expr.targ     <- get(target_data)$rankz
expr.med      <- get(mediator_data)$rankz
  
  
  
  



# Get range of index
rng <- 1:nrow(mediation_data)
  
if(!is.null(chunk_number)){
   max_col = nrow(mediation_data)
   rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
   if(rng[length(rng)] > max_col) {
      
      rng = rng[1]:max_col
      
   }
}



mediation_data <- mediation_data[rng,] 
  







### Extract data
target.id         <- mediation_data$target.id
mediator.id       <- mediation_data$mediator.id
target.qtl.chr    <- mediation_data$target.qtl.chr
target.qtl.marker <- mediation_data$marker.id



### Allocate vectors to store info
n <- nrow(mediation_data)
causal.p       <- numeric(length = n)
reactive.p     <- numeric(length = n)
independent.p  <- numeric(length = n)
undecided.p    <- numeric(length = n)
best.model     <- character(length = n)
best.model.p   <- numeric(length = n)
best.triad     <- character(length = n)
best.triad.p   <- numeric(length = n)
best.mediation <- numeric(length = n)
best.LRmed     <- numeric(length = n)
best.LR        <- numeric(length = n)
t.d_t          <- numeric(length = n) 
m.d_m          <- numeric(length = n)
t.m_t          <- numeric(length = n) 
m.t_m          <- numeric(length = n) 
t.md_t.m       <- numeric(length = n)
t.md_t         <- numeric(length = n)









### Mediation test begin
for(i in 1:n){
    
    
    
    # Mediation Test
    med <- mediation_test(target = expr.targ[, target.id[i], drop = FALSE],
                          mediator = expr.med[,mediator.id[i], drop = FALSE],
                          annotation = annot.id.med[mediator.id[i], ],
                          driver = genoprobs[[target.qtl.chr[i]]][,,target.qtl.marker[i]],
                          kinship = K[[target.qtl.chr[i]]],
                          covar_tar = covar.targ,
                          covar_med = covar.med,
                          verbose = FALSE)
   
    med.best <- as.data.frame(med$test)
   
   
   
    ### Save results for ith index
    causal.p[i]       <- med$test[med$test$model %in% 'causal', 'pvalue', drop = TRUE]
    reactive.p[i]     <- med$test[med$test$model %in% 'reactive', 'pvalue', drop = TRUE]
    independent.p[i]  <- med$test[med$test$model %in% 'independent', 'pvalue', drop = TRUE]
    undecided.p[i]    <- med$test[med$test$model %in% 'undecided', 'pvalue', drop = TRUE]
    best.model[i]     <- as.character(med.best[which.min(med.best$pvalue),'model'])
    best.model.p[i]   <- med.best[which.min(med.best$pvalue),'pvalue']
    best.triad[i]     <- as.character(med$best$triad)
    best.triad.p[i]   <- med$best$pvalue
    best.mediation[i] <- med$best$mediation
    best.LRmed[i]     <- med$best$LRmed
    best.LR[i]        <- med$best$LR
    t.d_t[i]          <- med$fitsLR$t.d_t
    m.d_m[i]          <- med$fitsLR$m.d_m
    t.m_t[i]          <- med$fitsLR$t.m_t
    m.t_m[i]          <- med$fitsLR$m.t_m
    t.md_t.m[i]       <- med$fitsLR$t.md_t.m
    t.md_t[i]         <- med$fitsLR$t.md_t
   
   
    print(paste0(i,'out of',n))
}
  






### Save results
full_result <- cbind(mediation_data, data.frame(causal.p, reactive.p, independent.p, undecided.p,
                                                best.model, best.model.p, best.triad, best.triad.p,
                                                best.mediation, best.LRmed, best.LR,
                                                t.d_t, m.d_m, t.m_t, m.t_m,t.md_t.m,t.md_t))








### Save data as .rds
file_name <- c(target_id, mediator_id)      # Just for naming conventions
file_name <- gsub('_id','',file_name)
file_name <- gsub('gene','mrna', file_name)

if(!is.null(chunk_number)){
   saveRDS(full_result, paste0('attie_islet_',file_name[1],'_full_mediation_',file_name[2],'_chunk_',chunk_number,'.rds'))
}else{
   saveRDS(full_result, paste0('attie_islet_',file_name[1],'_full_mediation_',file_name[2],'.rds'))
}
