### Downlaod the individual TCGA expression for GBM and LGG and Related metadata
setwd("./GliomaGEM/")

library(biomaRt)
library(plyr)
library(TCGAbiolinks) 
library(SummarizedExperiment)
library(org.Hs.eg.db) 
library(DT)
library(data.table)
library(tidyverse)
FPKMtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}
# Downloading FPKM of GBM
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",#*Gene expression",
                  data.type = 'Gene Expression Quantification', #"Gene expression quantification",
                  #platform = "Illumina HiSeq",
                  #file.type  = 'HTSeq-FPKM',#"results",
                  workflow.type = "STAR - Counts",# "HTSeq - FPKM-UQ", 
                  #experimental.strategy = "RNA-Seq",
                  legacy = FALSE
                  )
# query <- GDCquery(
#   project = "TCGA-GBM",
#   data.category = "Gene expression",
#   data.type = "Gene expression quantification",
#   platform = "Illumina HiSeq", 
#   file.type  = "normalized_results",
#   experimental.strategy = "RNA-Seq",
#   legacy = TRUE
# )
GDCdownload(query, method = "api")
tcga_gbm_exp <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")
tcga_gbm_meta <- as.data.frame(colData(tcga_gbm_exp)) ## Downloading metadata
tcga_gbm_exp=assay(tcga_gbm_exp)
tcga_gbm_exp = as.data.frame(tcga_gbm_exp)
#saveing GBM FPKM
fwrite(tcga_gbm_exp,'./data/TCGBiolinks/TCGA_GBM_FPKM_TCGBiolinks.csv',row.names = TRUE)

# Convert FPKM to TPM
tcga_gbm_tpm <- apply(as.matrix(tcga_gbm_exp), 2, function(x) FPKMtoTPM(x))
tcga_gbm_tpm = as.data.frame(tcga_gbm_tpm)
fwrite(tcga_gbm_tpm,'./data/TCGBiolinks/TCGA_GBM_TPM_TCGBiolinks.csv',row.names = TRUE)


# Downloading FPKM of LGG
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",#*Gene expression",
                  data.type = 'Gene Expression Quantification', #"Gene expression quantification",
                  #platform = "Illumina HiSeq",
                  #file.type  = 'HTSeq-FPKM',#"results",
                  workflow.type = "STAR - Counts",#"HTSeq - FPKM",
                  #experimental.strategy = "RNA-Seq",
                  legacy = FALSE)

# query <- GDCquery(
#   project = "TCGA-LGG",
#   data.category = "Gene expression",
#   data.type = "Gene expression quantification",
#   platform = "Illumina HiSeq", 
#   file.type  = "normalized_results",
#   experimental.strategy = "RNA-Seq",
#   legacy = TRUE
# )

GDCdownload(query, method = "api")
tcga_lgg_exp <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")
tcga_lgg_meta <- as.data.frame(colData(tcga_lgg_exp)) ## Downloading metadata

tcga_lgg_exp=assay(tcga_lgg_exp)
#tcga_lgg_exp = as.data.frame(tcga_lgg_exp)
#fwrite(as.data.frame(tcga_lgg_exp),'./data/TCGA_LGG_FPKM_TCGBiolinks.csv',row.names = TRUE)

# Convert FPKM to TPM
tcga_lgg_tpm <- apply(as.matrix(tcga_lgg_exp), 2, function(x) FPKMtoTPM(x))
tcga_lgg_tpm <- data.frame(tcga_lgg_tpm)
#fwrite(tcga_lgg_tpm,'./data/TCGA_GBM_TPM_TCGBiolinks.csv',row.names = TRUE)

##### Merging TCGABiolinks metadata

tcga_gbm_meta <- tcga_gbm_meta[,!colnames(tcga_gbm_meta) %in% 'treatments']
tcga_lgg_meta <- tcga_lgg_meta[,!colnames(tcga_lgg_meta) %in% 'treatments']

#find intersected metadata
mergeCols <- intersect(colnames(tcga_lgg_meta),colnames(tcga_gbm_meta))
tcga_biolinks_meta<- rbind(tcga_gbm_meta[,mergeCols],tcga_lgg_meta[,mergeCols])

# Selecting the most relevant metadata
# others can be added
tcga_meta <- tcga_biolinks_meta[,c('barcode','sample','definition','paper_IDH.status','paper_X1p.19q.codeletion','primary_diagnosis',
                                   "paper_Histology","paper_Grade","paper_MGMT.promoter.status",
                                   "paper_Age..years.at.diagnosis." , "paper_Survival..months." ,"paper_Vital.status..1.dead.",
                                   "paper_Gender")]
tcga_meta[tcga_meta$definition=='Solid Tissue Normal','primary_diagnosis'] <- 'Solid Tissue Normal'
# Replacing NA
tcga_meta$paper_IDH.status <- addNA(tcga_meta$paper_IDH.status)
tcga_meta$paper_X1p.19q.codeletion <- addNA(tcga_meta$paper_X1p.19q.codeletion)
tcga_meta$primary_diagnosis <- addNA(tcga_meta$primary_diagnosis)
tcga_meta$primary_diagnosis = gsub('\\, ','_',tcga_meta$primary_diagnosis)
tcga_meta$primary_diagnosis = gsub('\\ ','_',tcga_meta$primary_diagnosis)

tcga_meta$paper_X1p.19q.codeletion = gsub('\\-','_',tcga_meta$paper_X1p.19q.codeletion)

### Stratifying the samples to 3 glioma subtypes according to WHO 2021 classification:
## Astrocytoma, IDH-mutant
## Oligodendroglioma, IDH-mutant and 1p/19q-codeleted
## Glioblastoma, IDH-wildtype

tcga_meta['WHO_2021'] ="NA_"
for(i in 1:nrow(tcga_meta)) {
  if (tcga_meta[i,"definition"]=='Solid Tissue Normal') { # Keeping the 5 control samples
    tcga_meta[i,'WHO_2021'] <- 'CTRL'
    # any sample with missing paper_Histology will be given "NA_" and can be removed
  } else if(is.na(tcga_meta[i,"paper_Histology"])  ) {
    tcga_meta[i,'WHO_2021'] <- 'NA_'
    # Mixed glioma (oligoastrocytoma) is added as a quality control , but it's considered poorly defined subtype
  } else if(tcga_meta[i,"paper_Histology"]== "oligoastrocytoma" ) {
    tcga_meta[i,'WHO_2021'] <- 'ODG_AST'
    
    # Since the main 3 glioma subytpes depends on IDH mutation, any sample with missing IDH status will be removed
  } else if(is.na(tcga_meta[i,"paper_IDH.status"])  ) {
    tcga_meta[i,'WHO_2021'] <- 'NA_'
    
    # Only IDH_wt will be selected for GBM
  } else if(tcga_meta[i,"paper_IDH.status"]=='WT'  && tcga_meta[i,"paper_Histology"]=='glioblastoma') {
    tcga_meta[i,'WHO_2021'] <- 'GBM_IDH_wt'
    #  IDH_mut GBM will be removed for simplification, but it may be added to 'AST_IDH_mut' as mentioned by WHO 2021
  } else if(tcga_meta[i,"paper_IDH.status"]=='Mutant'  && tcga_meta[i,"paper_Histology"]=='glioblastoma') {
    tcga_meta[i,'WHO_2021'] <- 'NA_'
    
  } else if (tcga_meta[i,"paper_IDH.status"]=='Mutant' && tcga_meta[i,"paper_Histology"]=='astrocytoma') {
    tcga_meta[i,'WHO_2021'] <- 'AST_IDH_mut'
    
  } else if (tcga_meta[i,"paper_Histology"]=='oligodendroglioma' && tcga_meta[i,"paper_IDH.status"]=='Mutant' && tcga_meta[i,"paper_X1p.19q.codeletion"]=='codel') {
    tcga_meta[i,'WHO_2021'] <- 'ODG_IDH_mut_Codel'  

  }
}

fwrite(tcga_biolinks_meta,'./data/TCGA_TCGBiolinks_metadata.csv',row.names = TRUE) # whole metadata
fwrite(tcga_meta,'./data/TCGA_TCGBiolinks_metadata_Summary.csv',row.names = TRUE) # Only selected features




tcga_meta <- read_csv('./data/TCGA_TCGBiolinks_metadata_Summary.csv')
# Number of samples in each subtype
table(tcga_meta$WHO_2021)

# Visualizing the metadata itself
summary(tcga_meta[,c('definition','paper_IDH.status','paper_X1p.19q.codeletion','primary_diagnosis')])
x <- table(tcga_meta[,c('definition','paper_IDH.status','paper_X1p.19q.codeletion','primary_diagnosis')], useNA = "always")
x<- as.data.frame(x)
x <- x[x$Freq>0,]

fwrite(x,'./data/TCGA_TCGBiolinks_metadata_statistics.csv',row.names = TRUE)

# Creating a Mosaoic plot of 3 features  definition + paper_IDH.status + primary_diagnosis
library(vcd)
library(MASS)
library(forcats)
library(viridis)

tcga_meta_ <- tcga_meta[,c('definition','paper_IDH.status','paper_X1p.19q.codeletion','primary_diagnosis')] #'definition'

tcga_meta_table <- structable(~definition+paper_IDH.status+primary_diagnosis,data=tcga_meta_ )

# Mosiac plot for definition + primary_diagnosis + paper_IDH.status
mosaic(~ definition + primary_diagnosis + paper_IDH.status,  data = tcga_meta_table,
       split_vertical = c(TRUE, FALSE, TRUE),
       labeling_args = list(rot_labels = c(bottom = 90, top = 90,left=0,right=90)),
       margins = c(left = 13, bottom = 5,top=10),
       offset_labels = c(left = 4.6,top=3.2,bottom=1),
       offset_varnames=c(left= 9))


library(ggmosaic)
flights <- fly  %>%
  filter(!is.na(do_you_recline), !is.na(rude_to_recline))

tcga_meta_ <- tcga_meta[,c('paper_Histology','paper_IDH.status','paper_X1p.19q.codeletion','WHO_2021')]
tcga_meta_ <- na.omit(tcga_meta_)
colnames(tcga_meta_) <- str_replace(colnames(tcga_meta_) ,'paper_','')

tcga_meta_$IDH.status[tcga_meta_$IDH.status=='WT'] <- 'IDH Wt '
tcga_meta_$IDH.status[tcga_meta_$IDH.status=='Mutant'] <- 'IDH Mut '

tcga_meta_$X1p.19q.codeletion[tcga_meta_$X1p.19q.codeletion=='non_codel'] <- 'Non Codel'
tcga_meta_$X1p.19q.codeletion[tcga_meta_$X1p.19q.codeletion=='codel'] <- 'Codel'

tcga_meta_$IDH.status <- factor(tcga_meta_$IDH.status,c('IDH Wt ','IDH Mut '))
tcga_meta_$X1p.19q.codeletion <- factor(tcga_meta_$X1p.19q.codeletion,c('Non Codel','Codel'))
tcga_meta_$Histology <- str_to_sentence(tcga_meta_$Histology )
tcga_meta_$Histology <- factor(tcga_meta_$Histology,rev(c('Astrocytoma','Glioblastoma','Oligodendroglioma','Oligoastrocytoma')))

tcga_meta_$WHO_2021

mosaic2_examp <- ggplot(data = tcga_meta_) +
  geom_mosaic(aes(x=product(IDH.status,Histology, X1p.19q.codeletion),
                  fill = Histology, alpha = IDH.status),offset = 0.02, na.rm = TRUE) + 
  geom_mosaic_text(aes(x = product(IDH.status,Histology, X1p.19q.codeletion)), na.rm = TRUE, repel = TRUE)+
  scale_fill_viridis(discrete = TRUE)+
  scale_alpha_manual(values =c(.6,.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size = 12),    
        axis.text.y = element_text(angle = 0,size = 12),
        axis.title.x = element_text(angle = 0,size = 14),
        axis.title.y = element_text(angle = 90,size = 14)) + 
  labs(y="Histology", x="IDH mutation status + 1p.19q.codeletion", 
       title = "Stratification of TCGA-GBM & TCGA-LGG patients\nbased on IDH mutation and 1p/19q codeletion")
mosaic2_examp

ggsave('Figure/Fig_XX_metadata.png', mosaic2_examp,units = 'in', width = 10,height = 8,dpi = 300)

tcga_meta_$IDH_Status_
tcga_meta_table <- structable(~paper_IDH.status+paper_Histology+paper_X1p.19q.codeletion,data=tcga_meta_ )
mosaic( paper_X1p.19q.codeletion~ paper_IDH.status+  paper_Histology,  data = tcga_meta_table,na.action = na.omit,
       split_vertical = c(TRUE, FALSE),color = TRUE,las = 1,
       labeling_args = list(rot_labels = c(bottom = 90, top = 0,left=0,right=90)),
       margins = c(left = 9, bottom = 5,top=4),
       offset_labels = c(left = 3,top=0,bottom=1),
       offset_varnames=c(left= 6,bottom=2.5))
library(ggpubr)
ggballoonplot(housetasks, fill = "value")+
  scale_fill_viridis_c(option = "C")

tcga_meta %>% group_by(paper_Histology,primary_diagnosis,paper_Grade) %>% 
  mutate(IDH_WT_Percentage = length(paper_IDH.status[paper_IDH.status=="WT"]) /length(paper_IDH.status[paper_IDH.status %in%c("Mutant","WT")]),
         Non_Codeletion_Percentage = length(paper_X1p.19q.codeletion[paper_X1p.19q.codeletion=="non_codel"]) /length(paper_X1p.19q.codeletion%in%c("non_codel","codel")),
         count = n()) %>%
  distinct(paper_Histology,primary_diagnosis,paper_Grade,count,IDH_WT_Percentage,Non_Codeletion_Percentage)-> tcga_meta_summ
tcga_meta_summ[,4:6] <- round(tcga_meta_summ[,4:6],2)
fwrite(tcga_meta_summ,'./data/TCGA_TCGBiolinks_metadata_statistics.csv',row.names = TRUE)

