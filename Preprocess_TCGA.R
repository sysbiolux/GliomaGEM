setwd("./GliomGEM")
library(DESeq2)
library(knitr)
library(arrayQualityMetrics)
library(devtools)

#Check if FPKM are well clustered in the PCA
# according to the 3 glioma types GBM, DDG, AST
tcga_meta = read.csv('data/TCGA/TCGA_mRNAseq_702_clinical.txt',sep='\t') 

#drop missing values from IDH values()
#tcga_meta = tcga_meta[complete.cases(tcga_meta[,c('IDH_mutation_status','X1p19q_codeletion_status')]), ]
# Make a new class column based on IDH and codeltion only
tcga_meta['TYPE'] =0
for(i in 1:nrow(tcga_meta)){
  
  if (is.na(tcga_meta[i,"IDH_mutation_status"]) ) {
    tcga_meta[i,'TYPE'] <- 'UNK'
  } else if (tcga_meta[i,"IDH_mutation_status"]=='Mutant' && is.na(tcga_meta[i,"X1p19q_codeletion_status"]))  {
    tcga_meta[i,'TYPE'] <- 'UNK'
  } else if(tcga_meta[i,"IDH_mutation_status"]=='Wildtype') {
    tcga_meta[i,'TYPE'] <- 'GBM'
  } else if (tcga_meta[i,"IDH_mutation_status"]=='Mutant' && tcga_meta[i,"X1p19q_codeletion_status"]=='Non-codel')  {
    tcga_meta[i,'TYPE'] <- 'AST'
  } else if (tcga_meta[i,"IDH_mutation_status"]=='Mutant' && tcga_meta[i,"X1p19q_codeletion_status"]=='Codel') {
    tcga_meta[i,'TYPE'] <- 'ODG'    
  }
}

library(tidyr)
library(ggplot2)
require(dplyr)
tcga_meta_cnt <- tcga_meta %>% count(Histology, TYPE)
# Group analysis of the assigned class and labeled Histology
png(filename="Figure/TCGA_Histology_Barplot.png", units="in", width=10, height=10, res=300)
ggplot(tcga_meta_cnt, aes(Histology, n, fill = TYPE)) + geom_col(position = "dodge")
dev.off()

#Read TCGA data
tcga = read.csv('data/TCGA/TCGA_mRNAseq_702.txt',sep='\t',row.names = 1) 


library("FactoMineR")
library("factoextra")
PC <- PCA(t(tcga), graph = FALSE)
fviz_pca_ind(PC, 
             col.ind = factor(tcga_meta$TYPE),
             #palette = pal_m,
             palette = "jco",
             #col.ind.sup =factor(tTCGA_meta$TYPE),
             pointshape = 20, pointsize = 2,
             geom =  c("point"), 
             #addEllipses = TRUE,
             #ellipse.type = "confidence",  ellipse.level=0.95,
             title='log2(normalized counts)')

P_X_PCA_1 <- fviz_pca_ind(PC, axes=c(1,2), col.ind = tcga_meta$TYPE, pointshape = 20, pointsize = 2,
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='Title')

png(filename="Figure/TCGA_PC1_PC2.png", units="in", width=10, height=10, res=300)
dev.off()

png(filename="Figs/TCGA_PC1_PC3.png", units="in", width=10, height=10, res=300)
fviz_pca_ind(PC, axes=c(1,3), col.ind = tcga_meta$TYPE, pointshape = 20, pointsize = 2, geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='Title')
dev.off()

png(filename="Figs/TCGA_PC2_PC3.png", units="in", width=10, height=10, res=300)
fviz_pca_ind(PC, axes=c(2,3), col.ind = tcga_meta$TYPE, pointshape = 20, pointsize = 2, geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='Title')
dev.off()
