library(dplyr)


### Read in the required data
load("~/Desktop/Attie/Islet Proteins/attie_islet_protein_284_qtl2_input.Rdata")
lod.peaks <- readRDS("~/Desktop/Attie/Islet Proteins/Additive/attie_islet_protein_284_rZ_lodpeaks_6.rds") # 8350 x 13
expr.annots <- readRDS("~/Desktop/Attie/Islet Proteins/annotated_uniprotID.rds")   # 8235 x 11
ensembl.version <- 91


### Changing covar.factors for QTL Viewer
covar.factors <- data.frame(column.names = covar.factors)  # 284 x 12
covar.factors$display.name <- gsub("([A-Z])([A-Z])", "\\1\\2 ", covar.factors$column.names)
covar.factors$display.name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", covar.factors$display.name, perl=TRUE)



### Removing rows in the Uniprot to Ensembl ID dataframe with NAs and Uniprot IDs mapping to multiple Ensembl gene IDs.
expr.annots <- expr.annots[complete.cases(expr.annots),]  # 7224 x 11
expr.annots <- expr.annots[nchar(expr.annots$gene_id) == 18, ] # 7160 x 11



### Making sure column names match in raw, norm, rankz data
stopifnot(colnames(norm) == colnames(raw))
stopifnot(colnames(norm) == colnames(rankz))



### Converting Uniprot IDs in the colum names to Ensembl protein IDs if there is one.
data_names <- colnames(norm)

for(i in 1:length(data_names)){
  
  name <- data_names[i]
  ensembl_id <- expr.annots[expr.annots$Majority.protein.IDs == name, 'protein_id']
  
  if(length(ensembl_id)!=0){
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
#raw <- raw[,grep('ENSMUSP', colnames(raw))]         # 284 x 5433 
norm <- norm[,grep('ENSMUSP', colnames(norm))]      # 284 x 5433 -> 284 x 4826
rankz <- rankz[,grep('ENSMUSP', colnames(rankz))]   # 284 x 5433 -> 284 x 4826



### Making sure column names match again
stopifnot(colnames(norm) == colnames(rankz))


### Converting the raw dataframe to matrix
raw <- as.matrix(raw)



### Chaning Uniprot IDs in the lod.peaks dataframe to Ensembl protein ID
for(i in 1:nrow(lod.peaks)){
  
  name <- lod.peaks$annot.id[i]
  ensembl_id <- expr.annots[expr.annots$Majority.protein.IDs %in% name, 'protein_id']
  
  if(length(ensembl_id) !=0){
    lod.peaks$annot.id[i] <- ensembl_id
  }else{
    lod.peaks$annot.id[i] <- lod.peaks$annot.id[i]
  }
  
}



### Removing rows in lod.peaks dataframe without an Ensembl protein ID
lod.peaks <- lod.peaks[grep('ENSMUSP', lod.peaks$annot.id),]     # 8216 x 13   ->  7267 x 13



### Making sure the start, end, middle, and strand are numerics
expr.annots$start <- as.numeric(expr.annots$start)
expr.annots$end <- as.numeric(expr.annots$end)
expr.annots$middle <- as.numeric(expr.annots$middle)
expr.annots$strand <- as.numeric(expr.annots$strand)



### Keeping rows in the Uniprot Annotation dataframe where the protein id is in the column names of raw, norm, and rankz
expr.annots <- expr.annots[expr.annots$protein_id %in% colnames(norm), ]   # 7160 x 11 -> 4826 x 11



### Making row names of annoted Uniprot ID to be the protein ID
rownames(expr.annots) <- expr.annots$protein_id


### Changing the original Uniprot ID column to a simpler name and rearranging the order of columns
colnames(expr.annots)[grep('Majority.protein.IDs',colnames(expr.annots))] <- 'uniprot_id'
expr.annots <- expr.annots[,c(3:ncol(expr.annots),1,2)]



### Add gene.symbol, gene.chr, gene.start, gene.end, and cis to lod.peaks data.frame
colnames(lod.peaks)[colnames(lod.peaks) %in% c('chr','pos')] <- c('qtl.chr','qtl.pos')
lod.peaks <- merge(lod.peaks, expr.annots[,c('protein_id','symbol','chr','start','end')], 
                   by.x = 'annot.id', by.y = 'protein_id', sort = FALSE)
colnames(lod.peaks)[colnames(lod.peaks) %in% c('symbol','chr','start','end')] <- c('gene.symbol','gene.chr','gene.start','gene.end')

lod.peaks = lod.peaks %>% mutate(cis = (gene.chr == qtl.chr) & (abs(gene.start - qtl.pos) <= 4))
lod.peaks <- lod.peaks[,c('annot.id','marker.id','lod','qtl.chr','qtl.pos','gene.chr','gene.start','gene.end', 'gene.symbol','cis','A','B','C','D','E','F','G','H')]


### Removing unnecessary data
rm(name, i, ensembl_id, data_names)


### Preparing dataset.protein
dataset.islet.proteins <- list(annots = expr.annots, 
                               covar = covar,
                               covar.factors = covar.factors,
                               datatype = datatype,
                               display.name = display.name,
                               lod.peaks = list(additive = lod.peaks),
                               norm = norm,
                               expr = rankz,
                               raw = raw,
                               samples = samples)

rm(expr.annots, covar, covar.factors, datatype, display.name, lod.peaks, norm, pheno.dict, rankz, raw, samples)
