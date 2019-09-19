### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)








### Load and get data
load('attie_islet_284_qtl_viewer_v2.RData')
prot_data <- 'dataset.islet.proteins.284'
mrna_data <- 'dataset.islet.mrna.284'

prot_annot <- get(prot_data)[['annot.protein']]
prot_covar <- get(prot_data)[['covar.matrix']]
prot_expr  <- get(prot_data)[['data']]$rz
mrna_annot <- get(mrna_data)[['annot.mrna']]
mrna_covar <- get(mrna_data)[['covar.matrix']]
mrna_expr  <- get(mrna_data)[['data']]$rz









### Find overlap
annots <- prot_annot %>% filter(gene.id %in% mrna_annot$gene.id)
prot_expr <- prot_expr[,annots$protein.id]
mrna_expr <- mrna_expr[,unique(annots$gene.id)]









### Help function
#     This function is used to regress out a covariate
#
#     Parameters:
#       expr  - expression matrix
#       covar - covariate matrix as created by model.matrix
#       out   - which covariate to regress out
regress_out <- function(expr, covar, out){
  
  result <- matrix(data = 0, nrow = nrow(expr), ncol = ncol(expr), dimnames = dimnames(expr))
  for(i in 1:ncol(expr)){
    
      temp <- cbind(data.frame(pheno = expr[,i]), covar)
    
      fit  <- lm(pheno ~., data = temp, na.action = na.exclude)
      coef <- coefficients(fit)[-c(which(is.na(coefficients(fit))), grep(out, names(coefficients(fit))))]
    
      stopifnot(which(is.na(residuals(fit))) == which(is.na(temp$pheno)))
      value <- coef[1] + (as.matrix(temp[,names(coef)[-1]]) %*% coef[-1]) + residuals(fit)  
    
      result[,i] <- value
  }
  
  result
}





### Regressing out sex before correlating
prot_expr2 <- regress_out(expr = prot_expr, covar = prot_covar, out = 'sex')
mrna_expr2 <- regress_out(expr = mrna_expr, covar = mrna_covar, out = 'sex')





### Correlation
pearson   <- apply(annots, 1, function(x) cor(prot_expr2[,x[3]], mrna_expr2[,x[4]], use = 'pairwise.complete.obs', method = 'pearson'))
pearson_p <- apply(annots, 1, function(x) cor.test(prot_expr2[,x[3]], mrna_expr2[,x[4]], use = 'pairwise.complete.obs', method = 'pearson')$p.value)
pearson_adj <- p.adjust(p = pearson_p, method = 'BH')









### Plot histogram
par(lwd = 4)
hist(x = pearson, breaks = seq(-0.4, 1, .2),
     border = 'darkseagreen3', col = 'grey97',
     angle = 45, density = 10, cex.axis = 1, cex.lab = 1.2,
     font.lab = 2,
     xlab = "Pearson's coefficent",
     main = 'Protein - Transcript Correlations')
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)


text(x = .5, y = 2000, 
     labels = bquote(mu ~ ""%~~%.(signif(mean(pearson), 2))), 
     cex = 1.5, col = 'black')



