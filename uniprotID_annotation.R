### This script takes the Attie Uniprot ID and ensembl transcript ID file
#       and converts them to gene_id/protein_id. Additional columns such as 
#       strand, chr, start, end, middle, symbol, nearest.marker.id are also included.
#
### Input: Attie's 'UniprotID_to_ENSMBL.txt' file.
#          A marker file. I used "marker_grid_0.02cM_plus.rds"
#
### Output: RData file containing the information listed above for each Uniprot and transcript id.
#
### Author: Duy Pham
### Date:   July 5, 2018
####################################################################################################


### Install required library packages
# source("https://bioconductor.org/biocLite.R")
# biocLite(c('AnnotationHub','rtracklayer'))
# install.packages(c('tidyr','dplyr','data.table))


### Load libraries
library(AnnotationHub)
library(rtracklayer)
library(tidyr)
library(dplyr)
library(data.table)


### Options
options(stringsAsFactors = FALSE)


### Set a directory to save the data (optional)
setwd("~/Desktop/pQTL_Project/data/")


### Read in the Uniprot ID and markers file
uniprot_file <- read.table('~/Desktop/pQTL_Project/data/UniprotID_to_ENSMBL.txt', sep = '\t', na.strings = 'N/A', header = TRUE)
markers <- readRDS("~/Desktop/pQTL_Project/data/islet_proteins_qtl/marker_grid_0.02cM_plus.rds")



### Line 47: Remove the last ';' at the end of each ENSEMBL.ID
### Line 49: Create a new row for each transcript ID in the ENSEMBL.ID column
### Line 51: Remove brackets and words in between the brackets if it is next to the transcript ID
### Line 54: Remove trailing white spaces
uniprot_file$ENSEMBL.ID <- substr(uniprot_file$ENSEMBL.ID, 1, nchar(uniprot_file$ENSEMBL.ID)-1)

uniprot_file <- uniprot_file %>% separate_rows('ENSEMBL.ID', sep = ';')

uniprot_file$ENSEMBL.ID <- sapply(uniprot_file$ENSEMBL.ID, 
                                        FUN = function(x) gsub("\\[[^]]*]", "", x))

uniprot_file$ENSEMBL.ID <- sapply(uniprot_file$ENSEMBL.ID, 
                                        FUN = function(x) trimws(x, which = 'right'))



### Load in the AnnotationHub database (db) for ENSEMBL IDs and info
hub <- query(AnnotationHub(), c('ensembl','gtf','mus musculus'))
hub <- hub[grep('GRCm38', hub$title)]
anno_hub <- hub[["AH60127"]]
anno_hub <- as.data.frame(anno_hub)

### Subsetting the AnnotationHub db based on 'type' to extract info
ensembl_t <- anno_hub[anno_hub$type == 'transcript',]
ensembl_g <- anno_hub[anno_hub$type == 'gene',]
ensembl_cds <- anno_hub[anno_hub$type == 'CDS',]



### Here I merge the gene_id column in the ensembl_t dataframe based on the transcript_id 
#     in both the ensembl_t and uniprot_file_split dataframe.
annots <- merge(uniprot_file, ensembl_t[,c('transcript_id','gene_id')], 
                  by.x = 'ENSEMBL.ID', by.y = 'transcript_id', all.x = TRUE, sort = FALSE)


### Next I merge the chromsome (seqnames), start, end, strand, and gene_name columns 
#     based on the gene_id in both the ensembl_g and anno_uniprot dataframe.
annots <- merge(annots, ensembl_g[,c('gene_id','seqnames','start','end','strand','gene_name')], 
                  by = 'gene_id', all.x = TRUE, sort = FALSE)


### Finally, I merge the protein_id based on the transcript_id 
#     in both the ensembl_cds and anno_uniprot dataframe.
ensembl_cds <- ensembl_cds[,c('transcript_id','protein_id')]
ensembl_cds <- unique(ensembl_cds)
annots <- merge(annots, ensembl_cds, by.x = 'ENSEMBL.ID', 
                  by.y = 'transcript_id', all.x = TRUE, sort = FALSE)



### Here I am replacing all ';' and '-' with '_' and '.', respectively in the Uniprot ID
#     to be consistent with the names in the QTL data. (97-98).
#
### Next, I am changing some of the columns as required by QTLViewer (99-104).
annots$Majority.protein.IDs <- gsub(';','_', annots$Majority.protein.IDs)  # Replace ';' with '_'
annots$Majority.protein.IDs <- gsub('-','.', annots$Majority.protein.IDs)  # Replace '-' with '.'
annots$strand <- as.character(annots$strand)  # Convert 'strand' column to character (was factor before)
annots[annots == "+"] <- 1                    # Convert "+" strand to 1 
annots[annots == "-"] <- -1                   # Convert "-" strand to 1 
annots$start <- annots$start / 1000000        # Convert start position to mbs (divide by 1000000)
annots$end <- annots$end / 1000000            # Convert end position to mbs (divide by 1000000)
colnames(annots)[c(1,4,8)] <- c('transcript_id','chr','symbol')   # Changing column name 'ENSEMBL.ID' to 'transcript_id'
                                                                  # Changing column name 'seqnames' to 'chr'
                                                                  # Changing column name 'gene_name' to 'symbol'



### Next I add a 'middle' column, which is the middle point between the 'start' and 'end' column.
annots$middle <- (annots$start + annots$end) / 2



### Here I extract the nearest marker to the Uniprot ID based on the 
#       chromosome (anno_uniprot) and 
#       minimum distance between the pos (markers) and middle point (anno_uniprot).
for(i in 1:nrow(annots)){
  
  ### Subset markers dataframe based on the chromosome of the i-th row in the anno_uniprot dataframe.
  submarkers <- subset(markers, chr == annots[i,'chr'])
  
  ### If there are no markers, store NA to the i-th row in the 'neareast.marker.id' column
  #     else, store the marker ID with the minimum distance to the middle point.
  
  if(nrow(submarkers) == 0){
    
    annots[i,'nearest.marker.id'] <- NA   
    
  }else{
    
    min_marker_index <- which.min(abs(annots[i,'start'] - submarkers$pos))
    annots[i,'nearest.marker.id'] <- submarkers[min_marker_index,'marker']
    
  }
}


                                  
### Here I merge all rows with the same Uniprot IDs as one. If a Uniprot ID has multiple rows, 
#     the value in each row for a particular column will be separated by ';'.
annots <- annots %>% group_by(Majority.protein.IDs) %>%
                     summarise(transcript_id = paste0(unique(transcript_id), collapse=";"),
                                protein_id = paste0(unique(protein_id), collapse=";"),
                                gene_id = paste0(unique(gene_id), collapse=";"),
                                symbol = paste0(unique(symbol), collapse =";"),
                                chr = paste0(unique(chr), collapse=";"),
                                start = paste0(unique(start), collapse=";"),
                                end = paste0(unique(end), collapse=";"),
                                strand = paste0(unique(strand), collapse=";"),
                                middle = paste0(unique(middle), collapse=";"),
                                nearest.marker.id = paste0(unique(nearest.marker.id), collapse=";"))



### Replacing character "NA" with NA and keeping the information before ;NA
#       Example: "protein_id1;NA" -> protein_id1
annots <- as.data.frame(annots)
annots[annots == "NA"] <- NA
annots <- as.data.frame(apply(annots, 2, FUN = function(x) gsub(';NA','',as.character(x))))
        
                              
                              
### Removing unnecessary data
rm(uniprot_file, ensembl_cds,ensembl_g,ensembl_t,hub,anno_hub, markers, submarkers, i, min_marker_index)
                              
                              
                              
### Save the data to .rds file                         
saveRDS(annots, 'annotated_uniprotID.rds')


