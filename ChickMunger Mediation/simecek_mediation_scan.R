######################################################################################################################
#
#   This script performs mediation scans given two QTL Viewer formatted datasets using
#        Petr Simecek's intermediate package: install_github('simecek/intermediate')
#
#
#   Input:
#      1.) Load in the QTL Viewer environment
#      2.) target_id    :  Name of the column in annots that will be used to select the targets
#      3.) target_data  :  dataset.* list where annots contain the column name in [2]
#      4.) mediator_id  :  Name of the column in annots that will be used to select the mediators
#      5.) mediator_data:  dataset.* list where annots contain the column name in [4]
#      6.) chunk_number :  numeric value indicating which portion of the lod peaks table to break. Need a fixed chunk_size
#      7.) chunk_size   :  numeric value indicating the size to break the lod peaks table
#
#
#
#   Output:
#      1.) Dataframe containing the target's info and each of the mediator's info as a one long string separated by ','
#
#
#
#   Author: Duy Pham
#   Date:   July 10, 2018
#   E-mail: duy.pham@jax.org
#
######################################################################################################################

### Load required library packages
options(stringsAsFactors = FALSE)
library(intermediate)
library(dplyr)



### Command line arguments/variables to change
args <- commandArgs(trailingOnly = TRUE)
load(args[1])
target_id     <- args[2]
target_data   <- args[3]
mediator_id   <- args[4]
mediator_data <- args[5]
chunk_number  <- as.numeric(args[6])
chunk_size    <- as.numeric(args[7])





# Get data
annot.id.targ <- get(target_data)$annots                                                 # Get target's annotation dataframe in dataset.*        (target_data)
annot.id.targ <- annot.id.targ %>% dplyr::rename(pos = start) %>% as.data.frame()        # Rename target_id column to id for mediation_scan     
annot.id.med  <- get(mediator_data)$annots                                               # Get mediator's annotation dataframe in dataset.*      (mediator_data)
annot.id.med  <- annot.id.med %>% dplyr::rename(pos = start) %>% as.data.frame()         # Rename mediator_id column to id for mediation_scan
covar.targ    <- get(target_data)$covar                                                  # Get target's covariate matrix in dataset.*            (target_data)
covar.med     <- get(mediator_data)$covar                                                # Get mediator's covariate matrix in dataset.*          (mediator_data)
expr.targ     <- get(target_data)$rankz                                                  # Get target's rankZ matrix in dataset.*                (target_data)
expr.med      <- get(mediator_data)$rankz                                                # Get mediator's rankZ matrix in dataset.*              (mediator_data)
lod.peaks     <- get(target_data)$lod.peaks$additive                                     # Get target's LOD peaks table in dataset.*             (target_data)








# Get range of index if chunk_numer and chunk_size is specified
rng <- 1:nrow(lod.peaks)

if(!is.null(chunk_number)){
  max_col = nrow(lod.peaks)
  rng = ((chunk_number - 1) * chunk_size + 1):(chunk_number * chunk_size)
  
  if(rng[length(rng)] > max_col) {
    rng = rng[1]:max_col
  }
}



lod.peaks <- lod.peaks[rng,]




# Get the lodpeaks info of the target (used for the mediation_scan)
qtl.annot  <- lod.peaks$annot.id
qtl.chr    <- lod.peaks$qtl.chr
qtl.marker <- lod.peaks$marker.id









# Create vectors to store target info
target.qtl.pos <- lod.peaks$qtl.pos
target.lod     <- lod.peaks$lod
target.chr     <- lod.peaks$gene.chr
target.start   <- lod.peaks$gene.start
target.end     <- lod.peaks$gene.end
target.symbol  <- lod.peaks$gene.symbol
target.cis     <- lod.peaks$cis


# Create empty vectors to store mediator info and mediation lod score
n = nrow(lod.peaks)
mediator.id     <- character(length = n)
mediator.chr    <- character(length = n)
mediator.start  <- character(length = n)
mediator.end    <- character(length = n)
mediator.symbol <- character(length = n)
mediation.lod   <- character(length = n) 








for(i in 1:n){
  
    # Mediation: Y ~ Q  + M + covar
    med <- mediation.scan(target = expr.targ[, qtl.annot[i], drop = FALSE],
                          mediator = expr.med,
                          annotation = annot.id.med,
                          qtl.geno = genoprobs[[qtl.chr[i]]][,,qtl.marker[i]],
                          covar = covar.targ,
                          method = "double-lod-diff",
                          verbose = FALSE)
  
  
  
    med <- med[med[,mediator_id] != qtl.annot[i],]
  
  
    # Save results
    mediator.id[i]     <- paste0(med[,mediator_id], collapse = ',')
    mediator.chr[i]    <- paste0(med$chr, collapse = ',')
    mediator.start[i]  <- paste0(med$pos, collapse = ',')
    mediator.end[i]    <- paste0(med$end, collapse = ',')
    mediator.symbol[i] <- paste0(med$symbol, collapse = ',')
    mediation.lod[i]   <- paste0(med$LOD, collapse = ',')
  
  
  
    print(paste0(i,'out of', n))
}





mediation_results <- data.frame(target.id = qtl.annot,           # Save target id
                                target.qtl.chr = qtl.chr,        # Save target's qtl chromosome
                                target.qtl.pos,                  # Save target's qtl position within chromosome
                                marker.id = qtl.marker,          # Save the name of the qtl marker
                                target.lod,                      # Save the LOD score of the target at the qtl marker
                                target.chr,                      # Save the target's chromosome location
                                target.start,                    # Save the target's start position
                                target.end,                      # Save the target's end position
                                target.symbol,                   # Save the target's symbol
                                target.cis,                      # Save the logical value (TRUE/FALSE) indicating whether the target's qtl is cis or distal
                                mediator.id,                     # Save the name of all mediator as one string separated by ','
                                mediator.chr,                    # Save the chromosome of all mediator as one string separated by ','
                                mediator.start,                  # Save the start position of all mediator as one string separated by ','
                                mediator.end,                    # Save the end position of all mediator as one string separated by ','
                                mediator.symbol,                 # Save the symbol of all mediator as one string separated by ','
                                mediation.lod)                   # Save the mediation score from all mediator as one string separated by ','








### Save data as .rds
file_name <- c(target_id, mediator_id)      # Just for naming conventions
file_name <- gsub('_id','',file_name)
file_name <- gsub('gene','mrna', file_name)

if(!is.null(chunk_number)){
  saveRDS(mediation_results, paste0('attie_islet_',file_name[1],'_simecek_intermediate_',file_name[2],'_chunk_',chunk_number,'.rds'))
}else{
  saveRDS(mediation_results, paste0('attie_islet_',file_name[1],'_simecek_intermediate_',file_name[2],'.rds'))
}


