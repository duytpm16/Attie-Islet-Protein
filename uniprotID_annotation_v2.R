options(stringsAsFactors = FALSE)
library(ensimplR)
library(xlsx)
library(dplyr)
setwd('~/Desktop')




uniprot <- read.delim("~/Desktop/Attie Mass Spectrometry/Islet/Proteins/UniprotID_to_ENSMBL.txt", na.strings = 'N/A')
ensembl <- readRDS("~/Desktop/Ensembl/Ensembl 94/reference_protein_info_ensembl_94.rds")




### Writing all uniprot ID to vector and then to .xlsx file
all_uniprot <- vector()
for(i in 1:nrow(uniprot)){

    all_uniprot <- c(all_uniprot, strsplit(uniprot$Majority.protein.IDs[i], split = ';')[[1]])

}
write.xlsx(all_uniprot, 'all_unprot.xlsx')





### Data after running uniprot conversion tool
uniprot_to_proteinID <- read.table('https://www.uniprot.org/mapping/M20190207216DA2B77BFBD2E6699CA9B6D1C41EB2049ECCO.tab', header = TRUE)
ensembl$protein_id    <- gsub("\\..*","",ensembl$protein_id)
ensembl$gene_id       <- gsub("\\..*","",ensembl$gene_id)
ensembl$transcript_id <- gsub("\\..*","",ensembl$transcript_id)
uniprot_to_proteinID <- merge(uniprot_to_proteinID, ensembl, by.x = 'To', by.y = 'protein_id', all.x = TRUE)


rownames(uniprot)     <- uniprot$Majority.protein.IDs
uniprot$protein_id    <- ''
uniprot$gene_id       <- ''
uniprot$transcript_id <- ''
for(i in 1:nrow(uniprot_to_proteinID)){
    
    # Get index in uniprot dataframe where uniprot ID exists
    index <- grep(uniprot_to_proteinID$From[i], uniprot$Majority.protein.IDs)
    
    
    # Making sure to get the right location if it exists in multiple. Some Uniprot are like QGHE1Z for one row and QGHE1Z-2 in another. Grep will pick both of these indices up
    if(length(index) > 1){
       s <- c()
       for(j in index){
           if(uniprot_to_proteinID$From[i] %in% strsplit(uniprot$Majority.protein.IDs[j],split = ';')[[1]]){
              s <- c(s,j)
           }
       }
       index <- unique(s)
    }
    
    
    
    # Concatenating the protein_id, gene_id, and transcript_id to location of uniprot_id
    if(uniprot[index, 'protein_id'] == ''){
       uniprot[index, 'protein_id']    <- uniprot_to_proteinID$To[i]
       uniprot[index, 'gene_id']       <- uniprot_to_proteinID$gene_id[i]
       uniprot[index, 'transcript_id'] <- uniprot_to_proteinID$transcript_id[i]
    }else{
       uniprot[index, 'protein_id']    <- paste0(uniprot[index, 'protein_id'],    ';', uniprot_to_proteinID$To[i])
       uniprot[index, 'gene_id']       <- paste0(uniprot[index, 'gene_id'],       ';', uniprot_to_proteinID$gene_id[i])
       uniprot[index, 'transcript_id'] <- paste0(uniprot[index, 'transcript_id'], ';', uniprot_to_proteinID$transcript_id[i])
    }
  
}





### Get unique protein, gene, and transcript IDs
for(i in 1:nrow(uniprot)){
    
    uniprot$gene_id[i]       <- paste0(unique(strsplit(uniprot$gene_id[i],       split = ';')[[1]]), collapse = ';')
    uniprot$protein_id[i]    <- paste0(unique(strsplit(uniprot$protein_id[i],    split = ';')[[1]]), collapse = ';')
    uniprot$transcript_id[i] <- paste0(unique(strsplit(uniprot$transcript_id[i], split = ';')[[1]]), collapse = ';')

}




### Filter out gene IDs with NA annotations to try and recover some uniprot annotations
multi_geneID <- uniprot[uniprot$protein_id != '',]
multi_geneID <- multi_geneID[grepl(';', multi_geneID $gene_id),]
for(i in 1:nrow(multi_geneID)){
  
    # Get the gene ID info 
    temp <- batchGenes(as.list(strsplit(multi_geneID$gene_id[i], split =';')[[1]]), version = 94)
    
    
    # Find which has NA annotations
    na.i  <- which(is.na(temp$symbol))
    na.id <- temp$identifier[na.i]
    gm.i  <- grep('Gm', temp$symbol)
    gm.id <- temp$identifier[gm.i]

    
    # Filter out Gene ID with NA annotations. Most likely due to the annotation not being part of the autosome and sex chromosomes
    if(length(na.id) != 0){
       for(j in 1:length(na.id)){
           multi_geneID$gene_id[i] <- gsub(na.id[j], '', multi_geneID$gene_id[i], fixed = TRUE)
       }
       count <- strsplit(multi_geneID$gene_id[i], split = ';')[[1]]
       if(length(grep('ENSMUSG', count)) == 1){
          multi_geneID$gene_id[i] <- gsub(';', '', multi_geneID$gene_id[i])
       }
    }
    if(length(gm.id) != 0){
      for(j in 1:length(gm.id)){
        multi_geneID$gene_id[i] <- gsub(gm.id[j], '', multi_geneID$gene_id[i], fixed = TRUE)
      }
      count <- strsplit(multi_geneID$gene_id[i], split = ';')[[1]]
      if(length(grep('ENSMUSG', count)) == 1){
        multi_geneID$gene_id[i] <- gsub(';', '', multi_geneID$gene_id[i])
      }
    }
    
    
    if(length(unique(temp$symbol)) == 1){
       if(length(unique(temp$chromosome)) == 1){
          ver.i  <- which.min(as.numeric(temp$gene_id_version))
          ver.id <- temp$identifier[ver.i]
          multi_geneID$gene_id[i] <- gsub(ver.id[j], '', multi_geneID$gene_id[i], fixed = TRUE)
          count <- strsplit(multi_geneID$gene_id[i], split = ';')[[1]]
          if(length(grep('ENSMUSG', count)) == 1){
             multi_geneID$gene_id[i] <- gsub(';', '', multi_geneID$gene_id[i])
          }
       }
    }
    
    
}




new_singleID <- multi_geneID[!grepl(';', multi_geneID$gene_id),]
for(i in 1:nrow(new_singleID)){
    uniprot[new_singleID$Majority.protein.IDs[i],] <- new_singleID[i,]
}
















gene_info <- batchGenes(as.list(unique(uniprot$gene_id)), version = 94)
new_uniprot <- merge(uniprot, gene_info, by.x = 'gene_id', by.y = 'gene_id', all.x = TRUE)
new_uniprot <- new_uniprot %>%
                           dplyr::rename(chr = chromosome,
                                         uniprot_id = Majority.protein.IDs,
                                         description = name) %>%
                            dplyr::mutate(start = as.numeric(start),
                                          end = as.numeric(end),
                                          middle = (start + end) / 2,
                                          start = start / 1e6,
                                          end = end / 1e6,
                                          middle = middle / 1e6) %>%
                            dplyr::select(protein_id, transcript_id,gene_id, symbol,
                                          description, chr, start,middle, end, strand, synonyms,
                                          uniprot_id, entrez_id) %>%
                            arrange(protein_id)


saveRDS(new_uniprot, 'annotated_uniprots_v2.rds')



