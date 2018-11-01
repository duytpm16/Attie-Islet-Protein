options(stringsAsFactors = FALSE)
library(dplyr)


### Read in the required data
load("~/Desktop/Attie Final/Islet/Proteins/attie_islet_protein_qtl2_input.RData")
lod.peaks <- readRDS("~/Desktop/Attie Final/Islet/Proteins/Additive/attie_islet_protein_rZ_lodpeaks_6.rds") # 8856 x 13
expr.annots <- readRDS("~/Desktop/Attie Final/Islet/Proteins/annotated_uniprotID.rds")                      # 8235 x 11
ensembl.version <- 91


### Changing covar.factors for QTL Viewer
covar.factors <- data.frame(column.names = covar.factors)   # 373 x 15
covar.factors$display.name <- gsub("([A-Z])([A-Z])", "\\1\\2 ", covar.factors$column.names)
covar.factors$display.name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", covar.factors$display.name, perl=TRUE)



### Removing rows in the Uniprot to Ensembl ID dataframe with NAs and Uniprot IDs mapping to multiple Ensembl gene IDs.
expr.annots <- expr.annots[complete.cases(expr.annots),]         # 7224 x 11
expr.annots <- expr.annots[-grep(';',expr.annots$gene_id), ]     # 7160 x 11



### Making sure column names match in raw, norm, rankz data
stopifnot(colnames(norm) == colnames(raw))
stopifnot(colnames(norm) == colnames(rankz))



### Converting Uniprot IDs in the colum names to Ensembl protein IDs if there is one.
data_names <- colnames(norm)

for(i in 1:length(data_names)){
  
  ensembl_id <- expr.annots[expr.annots$Majority.protein.IDs %in% data_names[i], 'protein_id']
  
  if(length(ensembl_id) != 0){
     data_names[i] <- ensembl_id
  }else{
     data_names[i] <- data_names[i]
  }
}

colnames(raw) <- data_names
colnames(norm) <- data_names
colnames(rankz) <- data_names



### Making sure column names match again
stopifnot(colnames(norm) == colnames(raw))
stopifnot(colnames(norm) == colnames(rankz))



### Removing columns without Ensembl protein IDs
norm <- norm[,grep('ENSMUSP', colnames(norm))]       # 373 x 5433 -> 373 x 4808
rankz <- rankz[,grep('ENSMUSP', colnames(rankz))]    # 373 x 5433 -> 373 x 4808



### Making sure column names match again
stopifnot(colnames(norm) == colnames(rankz))


### Converting the raw dataframe to matrix
raw <- as.matrix(raw)



### Chaning Uniprot IDs in the lod.peaks dataframe to Ensembl protein ID
for(i in 1:nrow(lod.peaks)){
  
    name <- lod.peaks$annot.id[i]
    ensembl_id <- expr.annots[expr.annots$Majority.protein.IDs %in% name, 'protein_id']
  
    if(length(ensembl_id) != 0){
       lod.peaks$annot.id[i] <- ensembl_id
    }else{
       lod.peaks$annot.id[i] <- lod.peaks$annot.id[i]
    } 
}



### Removing rows in lod.peaks dataframe without an Ensembl protein ID
lod.peaks <- lod.peaks[grep('ENSMUSP', lod.peaks$annot.id),]     # 8856 x 13   ->  7857 x 13



### Making sure the start, end, middle, and strand are numerics
expr.annots$start  <- as.numeric(expr.annots$start)
expr.annots$end    <- as.numeric(expr.annots$end)
expr.annots$middle <- as.numeric(expr.annots$middle)
expr.annots$strand <- as.numeric(expr.annots$strand)



### Keeping rows in the Uniprot Annotation dataframe where the protein id is in the column names of raw, norm, and rankz
expr.annots <- expr.annots[expr.annots$protein_id %in% colnames(norm), ]   # 7160 x 11 -> 4808 x 11



### Making row names of annoted Uniprot ID to be the protein ID
rownames(expr.annots) <- expr.annots$protein_id



### Changing the original Uniprot ID column to a simpler name and rearranging the order of columns
expr.annots <- expr.annots %>% dplyr::rename(uniprot_id = Majority.protein.IDs) %>%
  select(protein_id, gene_id, symbol, chr, start, end, strand, middle, 
         nearest.marker.id, uniprot_id, transcript_id)




### Add gene.symbol, gene.chr, gene.start, gene.end, and cis to lod.peaks data.frame, then rearrange
lod.peaks <- lod.peaks %>% dplyr::rename(qtl.chr = chr, qtl.pos = pos)
lod.peaks <- merge(lod.peaks, expr.annots[,c('protein_id','symbol','chr','start','end')], 
                   by.x = 'annot.id', by.y = 'protein_id')
lod.peaks <- lod.peaks %>% dplyr::rename(gene.symbol = symbol, gene.chr = chr, gene.start = start, gene.end = end) %>%
  mutate(cis = (gene.chr == qtl.chr) & (abs(gene.start - qtl.pos) <= 4)) %>%
  select(annot.id,marker.id,lod,qtl.chr,qtl.pos,gene.chr,gene.start,gene.end,gene.symbol,cis,A,B,C,D,E,F,G,H)






### Preparing dataset.islet.proteins for QTL Viewer
dataset.islet.proteins <- list(annots = expr.annots, 
                               covar = covar,
                               covar.factors = covar.factors,
                               datatype = datatype,
                               display.name = display.name,
                               lod.peaks = list(additive = lod.peaks),
                               expr = norm,
                               rankz = rankz,
                               raw = raw,
                               samples = samples)



### Save data
rm(expr.annots, covar, covar.factors, datatype, display.name, lod.peaks, norm, 
   pheno.dict, rankz, raw, samples, name, i, ensembl_id, data_names)

save.image("~/Desktop/Attie Final/Islet/Proteins/attie_islet_protein_rZ_qtl_viewer.RData")
