library(qtl2)
library(qtl2convert)
library(dplyr)




### This function finds the nearest marker id to the 'chr' and 'start' position of the annotated protein/gene 
#      and returns annots with a new column 'nearest.marker.id'
#
find_nearest_marker <- function(id, annots, markers_df){
  
  
  # Creating new dataframe and find all unique chromosomes in markers
  temp_df <- data_frame()
  columns = c(id,'chr','start')
  
  
  # Subset annots and markers by chromosome and finds marker closest to the start position of the annotated id.
  for(i in unique(markers$chr)){
    
    sub_markers <- subset(markers_df, chr == i)
    
    sub_annot <- annots %>% 
      subset(chr == i) %>%
      group_by_at(vars(one_of(columns))) %>% 
      mutate(nearest.marker.id = sub_markers[which.min(abs(start - sub_markers$pos)),'marker'])
    
    temp_df <- bind_rows(temp_df, sub_annot)
  }
  
  
  
  # return annots with new column: 'nearest.marker.id'
  return(as.data.frame(temp_df))
  
  
  
}# find_nearest_marker





### This function runs a cis qtl scan for each annotated protein/gene. Returns a dataframe of the cis LOD scores if 'nearest.marker.id' column exists in annots, else
#      returns a list where the 1st element is a new annots dataframe containing a 'nearest.marker.id' column. The second element is a dataframe of the cis LOD scores.
#
cis_scan1 <- function(id, annots, genoprobs, pheno, kinship, addcovar, markers_df = NULL, intcovar = NULL){
       
      
      # If no 'nearest.marker.id' column exists, check to see if markers_df is supplied and run the find_nearest_marker function
      no.nearest.id <- FALSE
      
      if(is.null(annots$nearest.marker.id)){
         no.nearest.id <- TRUE 
        
         if(is.null(markers_df)){
            stop("Need markers dataframe since no columns with name, 'nearest.marker.id', exist in annots")
         
         }else{
            annots <- find_nearest_marker(id,annots,markers_df)
         }
      }
      
      
      
      annots <- annots[!is.na(annots$nearest.marker.id),]
      print0('Dimension of annotation dataframe: ', dim(annots))
      
      
      
      # Extracting some columns from annots so the code doesn't look messy
      chr = annots$chr
      nearest.marker <- annots$nearest.marker.id
      ids <- annots[,id, drop = TRUE]
      
      
      
      # Run QTL scan at nearest marker to protein/gene
      cis.lod <- matrix(NA, nrow = nrow(annots))
      for(i in 1:nrow(annots)){
          gp = genoprobs[,chr[i]]
          gp[[1]] = gp[[1]][,,nearest.marker[i], drop = FALSE]
        
          cis.lod[i,] <- scan1(genoprobs = gp, 
                                pheno = pheno[,ids[i], drop = FALSE],
                                kinship = K[[chr[i]]],
                                addcovar = addcovar,
                                intcovar = intcovar)
      }
      
      
      
      # Return data
      rownames(cis.lod) <- ids
      colnames(cis.lod) <- 'cis.lod'
      
      
      if(no.nearest.id){
        return(list(annots = as.data.frame(annots), cis.lod = cis.lod))
      
      }else{
        return(cis.lod)
      }
      
      
      
} #cis_scan1




### Example
mrna <- 'dataset.islet.mrna'
protein <- 'dataset.islet.proteins'

annot.id.protein <- get(protein)$annots
annot.id.mrna <- get(mrna)$annots


expr.protein <- get(protein)$expr
expr.mrna <- get(mrna)$expr

covar.protein <- get(protein)$covar
covar.mrna <- get(mrna)$covar

protein_cis_lod <- cis_scan1(id = 'protein_id', annots = annot.id.protein, genoprobs = genoprobs, pheno = expr.protein, kinship = K, addcovar = covar.protein, markers_df = markers)


mrna_cis_lod <- cis_scan1(id = 'gene_id', annots = annot.id.mrna, genoprobs = genoprobs, pheno = expr.mrna, kinship = K, addcovar = covar.mrna, markers_df = markers)



saveRDS(protein_cis_lod, 'attie_islet_protein_cis_lod_284.rds')
saveRDS(mrna_cis_lod, 'attie_islet_mrna_cis_lod_284.rds')
