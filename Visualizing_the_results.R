### Figure 1: Visualization of the FVA flux results of selected metabolites with literature evidence
### Figure 2: a) Predicted essential gene with common-essential and literature support
### Figure 2: b) Single drug targets
### Figure 3: Predicted combinations viability reduction in silico (a) and their targets (b)
### Figure 4: In vitro potency and CSF bio availability of the predicted drugs and combinations
### Figure 5: Xenografts data in screening databases (a) and literature (b)
### Figure 6: a) Summary of effective, ineffective, untested drugs in vitro, xenografts and clinical trials
### Figure 6: b) Survival measures of phase 2, 2-arms, clinical trials of the effective drugs 


setwd("./GliomGEM")
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}
wrap_text <- function(string, n) {
  spaces <- str_locate_all(string, " ")[[1]][,1]
  chars  <- nchar(string)
  for(i in 1:floor(chars/n)) {
    s <- spaces[which.min(abs(spaces - n*i))]
    substring(string, s, s) <- "\n "
  }
  return(string)
}
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
library(clipr)
#remotes::install_github("simmwill/coolors")
library(coolors)
library(data.table)
library(devtools)
library(DT)
#library(dplyr)
library(extrafont)
library(ggallin)
library(ggtext)
library(ggrepel)
library(ggplot2)
library(ggh4x)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplotify)
library(hrbrthemes)
library(knitr)
library(lemon)
#install.packages("pacman")
#install.packages("remotes")
#install.packages(c("here", "mice", "naniar"))
#remotes::install_github("emilelatour/lamisc")
library(lamisc) 
library(vroom)
library(magrittr)
library(openxlsx)
library(pheatmap)
library(patchwork)
library(readxl)
library(rrcov)
library(RColorBrewer)
library(stringi)
library(scales)
library(tidyr)
library(tidytext)
library(tidyverse)
library(tibble)
library(readr)
library(vegalite)
library(vroom)

set_breaks = function(limits) {
  seq(limits[1], limits[2], by = 1)
}

### Visualze TCGA metadata

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


### Figure 1:  Visualization of the FVA flux results and compare with literature evidence for sub type specific updates
fva <- read_csv('./Sample_models/Subtypes_Models_FVA.csv') # Fluxes of exchange reactions

fva %>% pivot_longer(cols = c( "minFlux_ast", "maxFlux_ast", "minFlux_gbm", "maxFlux_gbm" ,"minFlux_odg" ,"maxFlux_odg"),
                     names_to = 'Model',values_to = 'Flux') -> fva
fva %>% group_by(ex_rxns) %>%
  mutate(to_keep = if_else(is.na(Flux) | abs(Flux)<1,'remove','keep')) -> fva

#fva %>% group_by(ex_rxns) %>% filter(any(to_keep=='keep')) -> fva

#fva <- na.omit(fva)
#fva <- fva[fva$Flux>1e-4,]
unique(fva$ex_rxns)
fva$Subtype <- str_to_upper(str_split(fva$Model,'_',simplify = T)[,2])
fva$FluxType <- str_split(fva$Model,'Flux',simplify = T)[,1]
fva %>% group_by(ex_rxns) %>%
  mutate(average_flux = mean(abs(Flux), na.rm=TRUE)) -> fva
# Ceiling the measure to the 2nd decimal
fva$Flux_ <- lapply(as.numeric(fva$Flux) , function(x) signif.ceiling(x, 3))
fva$newname <- str_c(fva$Subtype,'_',fva$FluxType)
fva$FluxType[fva$FluxType=='min'] <- 'Min'
fva$FluxType[fva$FluxType=='max'] <- 'Max'
#fva$FluxType <- factor(fva$FluxType ,levels = c('Min','Max'))


met_list <- c("Glutamine","Thymidine","Glutamate",
              "4-Aminobutanoate","L-Phenylalanine","Myo-Inositol","Lactate")
fva$Subtype_rxn <- str_c(fva$Formula,fva$Subtype)
fva %>% group_by(Subtype_rxn) %>% mutate(ymin=min(Flux),ymax=max(Flux)) -> fva
fva$Subtype <- factor(fva$Subtype,c('GBM','AST','ODG'))
fva$Formula <- str_replace(fva$Formula,'  <=>','')
#fva %>% group_by(Formula, Subtype) %>%
#  mutate(Reaction=ifelse(sum(Flux)==0,'Inactive','Release')) ->fva
fva$Reaction <- NA
#fva$Reaction[fva$ymin=='NaN' & fva$ymax=='NaN' ]  <- 'Inactive'
#fva$Reaction[fva$ymin==0 & fva$ymax==0]  <- 'Inactive'

#fva$Reaction[fva$ymin=>-1e-12 & fva$ymax<=1e-12 ]  <- 'Inactive'

#fva$Reaction[fva$ymin=>0 & fva$ymax>0 ]  <-'Release' 3
#fva$Reaction[fva$ymax<=0 & fva$ymin<0 ]  <- 'Uptake'
#fva$Reaction[fva$ymin<0 & fva$ymax>0 ]  <-'Uptake/Release' 
#fva$Reaction[fva$ymin=>-1e-12 & fva$ymax<=1e-12 ]  <- 'Inactive'

fva$Reaction[fva$Flux==0 | fva$Flux =='NaN']  <- 'Inactive'
#fva$Reaction[fva$ymin=>0 & fva$ymax>0 ]  <-'Release' 3
fva$Reaction[fva$Flux <0 ]  <- 'Uptake'
fva$Reaction[fva$Flux >0 ]  <-'Release' 
fva$VMH_ID <- str_split(fva$ex_rxns,'EX_',simplify = T)[,2]
fva$VMH_ID <- str_split(fva$VMH_ID,'\\[',simplify = T)[,1]
fva$VMH_ID <- str_to_upper(fva$VMH_ID)

### Map VMH metaboolite IDs to super class
dico_met <- read_delim('./Generic_Models/metabolites_09_2021.csv',delim = '\t')
require(XML)
data <- xmlParse("csf_metabolites.xml")
#data <- xmlParse("hmdb_metabolites.xml")

xml_data <- xmlToList(data)
xml_data[5]$metabolite$taxonomy$direct_parent
length(xml_data$metabolite$accession)

class_df <- data.frame(metHMDBID=rep('',length(xml_data)), 
                       VMH_ID=NA, Name=NA,
                       super_class=NA, class=NA,sub_class=NA,molecular_framework=NA)

for(i in 1:length(xml_data)){
  class_df[i,'metHMDBID'] <-  as.character(xml_data[i]$metabolite$accession)
  class_df[i,'Name'] <-  xml_data[i]$metabolite$name
  
  vmh_id <- xml_data[i]$metabolite$vmh_id
  super_class <- xml_data[i]$metabolite$taxonomy$super_class
  class <- xml_data[i]$metabolite$taxonomy$class
  sub_class <- xml_data[i]$metabolite$taxonomy$sub_class
  molecular_framework <- xml_data[i]$metabolite$taxonomy$molecular_framework
  
  if (!is.null(vmh_id)){ class_df[i,'VMH_ID'] <-vmh_id }
  else{class_df[i,'VMH_ID'] <- NA    }

  if (!is.null(super_class)){ class_df[i,'super_class'] <-super_class }
  else{class_df[i,'super_class'] <- NA    }
  
  if (!is.null(class)){ class_df[i,'class'] <-class }
  else{class_df[i,'class'] <- NA    }
  
  if (!is.null(sub_class)){ class_df[i,'sub_class'] <-sub_class }
  else{class_df[i,'sub_class'] <- NA    }
  
  if (!is.null(molecular_framework)){ class_df[i,'molecular_framework'] <-molecular_framework }
  else{class_df[i,'molecular_framework'] <- NA    }
}
write_csv(class_df, 'Generic_Models/HMDB_Metabolie_Taxonomy.csv')

####

fva <- left_join(fva,class_df)
fva$super_class[is.na(fva$super_class)] <- 'Undefined'
unique(fva$Formula[is.na(fva$super_class)])

## Count the number of narriow-bound reactions
fva%>%  group_by(ex_rxns, Subtype)   %>% 
  mutate(rank_3_1 = Flux[FluxType=="Max"] -  Flux[FluxType=="Min"]) -> fva
n_distinct(fva$ex_rxns[fva$rank_3_1 >200])
min(fva$rank_3_1 )
max(na.omit(fva$rank_3_1 ))
min(na.omit(fva$rank_3_1 ))

n_distinct(fva$ex_rxns[fva$rank_3_1 >1e-4 & abs(fva$Flux) >250])
unique(fva$ex_rxns[fva$rank_3_1 >1e-4 & abs(fva$Flux) >250])

#dico_met <- left_join(dico_met,class_df)
# 
# dico_met$super_class <- NA
# dico_met$class <- NA
# dico_met$sub_class <- NA
# dico_met$molecular_framework<- NA
# for(i in 1:3){ #hmdb_ids
#   hmdb_id <- hmdb_ids[i]
#   x <- HmdbEntry(
#     prefix = "http://www.hmdb.ca/metabolites/",
#     id = hmdb_id,
#     keepFull = TRUE
#   )
#   y<- store(x)
#   super_class <- y$taxonomy$super_class
#   class<-  y$taxonomy$class
#   sub_class <- y$taxonomy$sub_class
#   molecular_framework <- y$taxonomy$molecular_framework
#   dico_met[dico_met$metHMDBID==hmdb_id, 'super_class'] <- super_class
#   dico_met[dico_met$metHMDBID==hmdb_id, 'class'] <- class
#   dico_met[dico_met$metHMDBID==hmdb_id, 'sub_class'] <- sub_class
#   dico_met[dico_met$metHMDBID==hmdb_id, 'molecular_framework'] <- molecular_framework
# }
#   
# data(hmdb1)
# y<- store(hmdb1)
# x <- diseases(hmdb1)
# names(store(hmdb1))
# biospecimens(hmdb1)
# tissues(hmdb1)
# hmdb1
# as.data.frame(hmdb1)
# x <- HmdbEntry(
#   prefix = "http://www.hmdb.ca/metabolites/",
#   id = "HMDB0000005",
#   keepFull = TRUE
# )
# biospecimens(x)
# y<- store(x)
# y$name
# y$taxonomy$super_class
# y$taxonomy$molecular_framework
# y$vmh_id

# Read the literature evidence of the exchange reaction to add in the plot
fva_literature <- read_excel('./Supplementary File 2.xlsx',sheet = "Table_S10_Exchange_Reactions")
colnames(fva_literature) <- c("Formula","Prediction","Literature","Subtype","Matching", "Cell Lines/ Tissue samples", "Reference"  )
#fva_literature$Literature <- wrap_text(fva_literature$Literature,n=25)
fva_literature$Formula <- str_split(fva_literature$Formula,'\ ',simplify = T)[,1]
fva_literature$FluxType <- c('Max','Max','Min','Max','Max','Max',"Max")

fva_selected <- fva[str_detect(fva$Formula,paste0(met_list,collapse = '|')),] 
fva_selected <- fva_selected[fva_selected$Formula!= "(R)-Lactate",]
fva_selected$Formula <- factor(fva_selected$Formula,fva_literature$Formula)
#fva_selected_2 <- fva[str_detect(fva$Formula,paste0(c("Glutamine","Thymidine"),collapse = '|')),] 
fva_selected <- left_join(fva_selected,fva_literature) # Joining, by = c("Formula", "Subtype")
fva_selected$Literature2 <- str_wrap(fva_selected$Literature , width = 40)
fva_selected$Subtype <- factor(fva_selected$Subtype,c("GBM","AST","ODG"))

fill = rev(c("cyan","cyan","cyan","cyan","cyan","cyan","magenta"))
muh_grob <- grid::rectGrob(
  x=0, y=1:7, gp=gpar(
    color='black', fill=fill, alpha=0.2))#rainbow(7)
fva_selected$Subtype_ <-str_c("i",fva_selected$Subtype)
fva_selected$Subtype_ <- factor(fva_selected$Subtype_,c("iGBM","iAST","iODG"))

P1 <- ggplot(fva_selected,aes(y=Formula,alpha=FluxType,x = abs(Flux)/10,group=Subtype_,shape=Reaction,
                        color=Subtype_,fill=Subtype_,size=FluxType))+
  geom_point(position = position_dodge(width = 0.5))+ #
  geom_linerange(aes(color = Subtype_,xmin=abs(ymin)/10,xmax=abs(ymax)/10), size = 1, alpha = 0.5,
                 position = position_dodge(width = 0.5))+
  scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'))+
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(4,6))+
  xlab("Absolute relative flux rate (%)") +
  ylab("Exchange reactions") +theme_classic()+
  theme(axis.text.x = element_text(angle = 0,size = 14),
        axis.text.y = element_text(angle = 0,size = 14,face='bold'),
        axis.title.y = element_text(size = 14),
        #panel.background = element_rect(fill = "grey92", colour = NA)
  ) + guides(size="none",alpha="none",fill=guide_legend(title = "Subtype\nmodel"),
             color=guide_legend(title = "Subtype\nmodel",legend.title=element_text(size=14),legend.label=element_text(size=14))) + 
  coord_cartesian(clip='off') +
  theme(axis.text.x = element_text(margin=margin(t=10),size =14),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12))+
  annotation_custom(grob=muh_grob, ymin = -0, ymax = 1, xmin = -24, xmax=8)
P1
P1_final <-P1 + xlim(0, 130) +
  geom_text_repel(aes(label=Literature2),size=3.8,alpha=0.9,
    #force        = 4,
    nudge_x      = 5,hjust=0,
    nudge_y      = 0.43,
    direction    = "both",
    fontface = 'bold',
    segment.colour = "black",
    box.padding = unit(0.35, "lines"),
    segment.size = 0.2
  )  

P1_final
ggsave(P1_final,filename="Figure/Figure_1_Exchange_Reactions.png", units="in", width=11, height=8.5)
#ggsave(P1_final,filename="Figure/Figure_1_Exchange_Reactions.pdf", units="in", width=11, height=8.5)

# 
# fva_literature$Literature2 <- sapply(strsplit(fva_literature$Literature, " "), function(x) {
#   spacePosition <- cumsum(nchar(x))
#   placeBreak <- spacePosition[which(diff(spacePosition %/% 30) == 1)] + 1
#   result <- paste(x, collapse = " ")
#   for(i in placeBreak) {
#     substring(result, i, i) <- "\n"
#   }
#   result
# })
# 
# ggplot(fva_literature,aes(y=Formula,x=10,label=Literature2))+
#   geom_label_repel()
# ## FVA for all reactions

fva %>% group_by(Formula) %>% 
  mutate(rank = sum(abs(Flux),na.rm = T))%>% 
  group_by(class) %>% 
  mutate(rank_2 = sum(abs(Flux),na.rm = T))-> fva

n_distinct(fva$class)
n_distinct(fva$sub_class)
unique(fva$molecular_framework)
CLASS <- sub_class
fva$Subtype_ <-str_c("i",fva$Subtype)
fva$Subtype_ <- factor(fva$Subtype_,c("iGBM","iAST","iODG"))

P_S1 <- ggplot(fva,aes(y=reorder_within(Formula,rank,rank_2),alpha=FluxType,x = abs(Flux)/10,
               group=Subtype_,
               shape=Reaction,
               color=class,fill=class,
              #color=sub_class,fill=sub_class,
              size=FluxType))+
  geom_point(position = position_dodge(width = 0.7))+ #
  geom_linerange(aes(xmin=abs(ymin)/10,xmax=abs(ymax)/10), size = 0.6, alpha = 0.5,
                 position = position_dodge(width = 0.7))+
  ggtitle("Uptake flux rate in the three glioma subtype models") +
  #scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'))+
  scale_y_reordered()+
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(2,3))+
  xlab("Absolute relative flux rate (%)") +
  ylab("Reactions") +theme_classic() +
  guides(size="none",alpha="none",
           color=guide_legend(legend.title=element_text(size=14),legend.label=element_text(size=14))) + 
  facet_grid(.~Subtype_,scales = "free")+ #Subtype
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 8),
  )
P_S1
ggsave(P_S1,filename="Figure/FigS_X_Subtype_Models_FVA.png", units="in", width=11, height=12)
ggsave(P_S1,filename="Figure/FigS_X_Subtype_Models_FVA.pdf", units="in", width=11, height=12)




### Figure 2: Essential gene with common-essential and literature support

### Add drug targets 
sko_df  <- read_csv('Integrated_Drug_Db/SKO_result_with_Effective_Targets.csv')
sko_df$TYPE2 <- "Single Drug Deletion"

dko_df  = read_csv('Integrated_Drug_Db/DKO_result_with_Effective_Targets.csv')
#colnames(dko_df) <- c('Drugs','Subtype','grRatio','DelRxns')
#dko_df$ <-'__'
dko_df$TYPE2 <- "Double Drug Deletion"
dko_df %>% separate_rows(Drugs,sep=' ;')%>%distinct() ->dko_df

sko_df <- rbind(sko_df,dko_df)
sko_df$TYPE2 <- factor(sko_df$TYPE2,c("Single Drug Deletion","Double Drug Deletion"))
sko_df <- sko_df[sko_df$Drugs !="fludarabine",]

unique(sko_df$Drugs)

## Convert ENTREZ to symbol sing in house dico disctionary
# sko_df_longer <- sko_df[,colnames(sko_df)!='All_Targets']
# cols <- c("Effective_Targets", "All_Targets")
# sko_df_longer %>% separate_rows(Effective_Targets,sep='; ')%>%distinct()%>%
#   mutate(entrez_id=str_split(Effective_Targets,"[.]",simplify = TRUE)[,1])->sko_df_longer

sko_df_longer <- sko_df[,colnames(sko_df)!='Effective_Targets']
cols <- c("Effective_Targets", "Effective_Targets")
sko_df_longer %>% separate_rows(All_Targets,sep='; ')%>%distinct()%>%
  mutate(entrez_id=str_split(All_Targets,"[.]",simplify = TRUE)[,1])->sko_df_longer

dico = read_csv('./Generic_Models/dico_201911.csv')
dico <- dico[,c('ENTREZ','SYMBOL')]
colnames(dico) <- c('entrez_id','SYMBOL')
dico <- na.omit(dico)
dico$entrez_id <- as.character(dico$entrez_id)
sko_df_longer <- left_join(sko_df_longer,dico)
sko_df_longer <- sko_df_longer[,c('Drugs','Subtype','TYPE2','SYMBOL','DelRxns')]

## Removing drug targets from  the individual drugs in the combination
integrated_db  <-  read_csv('Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv');
integrated_db <- integrated_db[integrated_db$name%in% sko_df_longer$Drugs,]
colnames(integrated_db)[1] <- 'Drugs'
sko_df_longer <- left_join(sko_df_longer,integrated_db[,c('Drugs','SYMBOL','database')],)

sko_df_longer <- distinct(sko_df_longer[!is.na(sko_df_longer$database),
                                        c('Drugs','Subtype','TYPE2','SYMBOL','DelRxns')])

#### Create a dataframe of gene, rxnname and rule
recon_pathways <- read_csv('Generic_Models/Recon3D_Rxn_Rules_Pathways.csv')
recon_pathways$SYMBOL <- recon_pathways$Rule
recon_pathways %>% separate_rows(SYMBOL,sep=' | ')%>%distinct()->recon_pathways_longer
recon_pathways_longer %>% separate_rows(SYMBOL,sep=' & ')%>%distinct()->recon_pathways_longer
recon_pathways_longer$SYMBOL <- str_replace_all(recon_pathways_longer$SYMBOL,'[(|)]','')
recon_pathways_longer <- recon_pathways_longer[!(recon_pathways_longer$SYMBOL %in% c("","&","|")),]
recon_pathways_longer <- recon_pathways_longer[! is.na(recon_pathways_longer$SYMBOL),]
recon_pathways_rule <- distinct(recon_pathways_longer[,c('SYMBOL',"Rule",'Pathway')])
recon_pathways_rxn <- distinct(recon_pathways_longer[,c("Rxn",'SYMBOL',"Rule",'Pathway')])

recon_pathways_longer <- distinct(recon_pathways_longer[,c('SYMBOL','Pathway')])

## Keep the targets of the deleted rxns
sko_df_longer_delrxn <- sko_df_longer
sko_df_longer_delrxn %>% separate_rows(DelRxns,sep='; ')%>%distinct()->sko_df_longer_delrxn
sko_df_longer_delrxn <- left_join(sko_df_longer_delrxn,recon_pathways_rxn,by=c('DelRxns'='Rxn')) #Rxn'='DelRxns',
sko_df_longer_delrxn %>% group_by(Drugs,Subtype,TYPE2) %>%
  summarise(SYMBOL = intersect(SYMBOL.x,SYMBOL.y)) -> sko_df_longer_delrxn
sko_df_longer <- sko_df_longer_delrxn
write_csv(sko_df_longer,"./Integrated_Drug_Db/Predicted_drug_targets_filtered.csv")

# ##  Removing targets with OR rule not covered
# sko_df_longer_rule <- left_join(recon_pathways_rule,sko_df_longer,by=c('SYMBOL'='SYMBOL')) #Rxn'='DelRxns',
# sko_df_longer_rule <- distinct(sko_df_longer_rule[! is.na(sko_df_longer_rule$TYPE2),])
# sko_df_longer_rule$OR_genes <- sko_df_longer_rule$Rule
# sko_df_longer_rule %>% separate_rows(OR_genes,sep=' | ')%>%distinct()->sko_df_longer_rule
# sko_df_longer_rule$OR_genes <- str_replace_all(sko_df_longer_rule$OR_genes,'[(|)]','')
# sko_df_longer_rule <- sko_df_longer_rule[!(sko_df_longer_rule$OR_genes %in% c("","&","|")),]
# 
# sko_df_longer_not_targeted <- sko_df_longer_rule %>% group_by(Drugs,Subtype,Rule) %>%
#   summarise(x = setdiff(OR_genes,SYMBOL),
#             Rule_targeted =   if_else(n_distinct(x)>0 ,'Yes','No')) 
#   
# #sko_df_longer_not_targeted <-  unique(sko_df_longer_rule[sko_df_longer_rule$OR_genes!= sko_df_longer_rule$SYMBOL,
# #                                                         c("Drugs","Subtype","Rule")])
# sko_df_longer_not_targeted <-  unique(sko_df_longer_not_targeted[, c("Drugs","Subtype","Rule")])
#                                                         
# sko_df_longer_not_targeted$Rule_targeted <- "No"
# sko_df_longer_targeted <- left_join(sko_df_longer_rule,sko_df_longer_not_targeted)
# #sko_df_longer_targeted <-  sko_df_longer_rule[!sko_df_longer_rule$Rule %in% sko_df_longer_not_targeted$Rule,]
# sko_df_longer_targeted <- sko_df_longer_targeted[is.na(sko_df_longer_targeted$Rule_targeted),]
# sko_df_longer <- distinct(sko_df_longer_targeted[,c('Drugs','Subtype','TYPE2','SYMBOL','DelRxns')])
# 
# ### Merge by rxn and gene 
# sko_df_longer %>% separate_rows(DelRxns,sep='; ')%>%distinct() -> sko_df_longer
# sko_df_longer
# sko_df_longer <- left_join(recon_pathways_longer,sko_df_longer,by=c('SYMBOL'='SYMBOL')) #Rxn'='DelRxns',
# sko_df_longer <- distinct(sko_df_longer[! is.na(sko_df_longer$TYPE2),])


####
genes_literature <- read_excel('./Supplementary File 2.xlsx',sheet = "Table_S11_Essential_Genes")

depmap_genes = read.csv("CRISPR/22Q1_Genes_Metadata.csv")

#data <- read.csv('./Sample_models/Essential_genes_UpSet_table_DMEM.csv',sep = ",") 
data <- read.csv('./Sample_models/Essential_genes_UpSet_table_.csv',sep = ",") 
depmap_scr_dep <- read_csv('./CRISPR/22Q1_DepMap_Chronos_Score_Probability.csv')

data
data%>%
  pivot_longer(!Row,names_to = "Genes", values_to = "Predicted_Essential") ->data_longer

data_longer$Row <- str_replace(data_longer$Row ,'Thiele2020_','')
data_longer$Row <- str_replace(data_longer$Row ,'__','_')
data_longer$Row <- str_replace(data_longer$Row,'_IDH_wt','')
data_longer$Row <- str_replace(data_longer$Row,'_IDH_mut','')
data_longer$Row <- str_replace(data_longer$Row,'_Codel','')
data_longer$Genes <- str_replace(data_longer$Genes,'X','')

data_longer%>%separate(Row, c("Model", "Data","Curation","Subtype"), "_") ->data_longer

# drop Human1 models
data_longer[str_detect(data_longer$Model,'Recon3D'),]-> data_longer
data_longer[str_detect(data_longer$Data,'Rahman2015'),]-> data_longer
data_longer[str_detect(data_longer$Curation,'CSF'),]-> data_longer
data_longer[!str_detect(data_longer$Subtype,'CTRL'),]-> data_longer

## number of common essentials in the predictions
colnames(data_longer) <- c("Model" ,   "Data"    , "Curation", "Subtype" , "Gene_ENTREZ" ,'Predicted_Essential')
depmap_genes <- depmap_genes[depmap_genes$Essentiality_Type!='Non-Essential',]
depmap_genes$Gene_ENTREZ <- as.character(depmap_genes$Gene_ENTREZ)
depmap_genes$Symbol <- str_replace(depmap_genes$Symbol," ","")
data_longer_ <- left_join(data_longer,depmap_genes)


## Cluster drugs by genes
sko_df_genes <- sko_df_longer %>%  group_by(SYMBOL,TYPE2) %>% 
  summarise(Drugs_ = paste0(unique(Drugs),collapse = "; " ))
sko_df_longer <- left_join(sko_df_longer, 
                           data_longer_[data_longer_$Predicted_Essential==1,c('Symbol',"Predicted_Essential")],
                           by=c('SYMBOL'='Symbol'))
sko_df_longer$Predicted_Essential[is.na(sko_df_longer$Predicted_Essential)] <- 0
sko_df_longer$SYMBOL_ <- sko_df_longer$SYMBOL
sko_df_longer$SYMBOL_[sko_df_longer$SYMBOL %in% c("CA7","CA6","CA5B","CA3","CA12","CA14","CA1","CA9","CA5A","CA4","CA2","CA8","CA13") ] <- "CA*"

sko_df_genes <- distinct(sko_df_longer[,c('SYMBOL','SYMBOL_','Drugs','Subtype','TYPE2',
                                          'Predicted_Essential')])#
sko_df_genes$Drugs <- firstup(sko_df_genes$Drugs)

sko_df_genes %>% group_by(SYMBOL_,TYPE2,Predicted_Essential,Subtype) %>%mutate(
  SYMBOL_id = str_length(paste0(unique(Drugs),collapse = "; " )),
  Drug_len = str_length(Drugs),n_drugs = n_distinct(Drugs),
  Drug_rank = ntile(rank(Drugs),n_distinct(Drugs)),
) ->sko_df_genes
#sko_df_genes$Subtype_ <- 'NA'
#sko_df_genes$Subtype_ <- sko_df_genes$Subtype
#sko_df_genes$Subtype_[sko_df_genes$Subtype_ %in% c('AST','ODG')] <- "AST/ODG"

sko_df_genes$Subtype <-  factor(sko_df_genes$Subtype,c('GBM','AST','ODG'))
#sko_df_genes$Subtype_ <-  factor(sko_df_genes$Subtype_,c('GBM',"AST/ODG"))

sko_df_genes$Predicted_Essential[sko_df_genes$Predicted_Essential==1] <- 'Essential drug targets'
sko_df_genes$Predicted_Essential[sko_df_genes$Predicted_Essential==0] <- 'Non-essential drug targets'
sko_placeholder__ <-  sko_df_genes  # A placeholder to add approved drug targets

sko_df_genes$Predicted_Essential <- factor(sko_df_genes$Predicted_Essential,c('Essential drug targets','Non-essential drug targets'))

#n_distinct(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion' & sko_df_genes$Predicted_Essential!= "Essential drug targets", "Drugs"])

sko_df_genes$Drugs_two_ch <- substr(sko_df_genes$Drugs ,1,2)
sko_df_genes$Drugs_one_ch <- substr(sko_df_genes$Drugs ,1,1)

mycols=c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gray",
                      "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
                      "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
                      "steelblue4", "darkturquoise", "green1", "yellow4",
                      "darkorange4", "yellow3", "brown", "skyblue3", "palegreen3", "khaki3", "#E31B4D",
                      "gold3","green2","orchid3","turquoise")

x <- unique(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',
                                  c("Drugs","Drugs_one_ch")])
# # Legend key width = maximum label * 1.10 to add some padding
# width <- unit(max(sapply(key_label, strwidth, units = "inches")) * 1.10, "in")
#unique(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion', "Pathway"])
#Y <- unique(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion', colnames(sko_df_genes)!="Pathway"])

#sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL %in% c("CA9","CA5A","CA4","CA2") ] <- "CA**"
#sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL %in% c("PI4KB","PI4KA","PI4K2A") ] <- "PI4K*"
#sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL %in% c("SLC7A1","SLC7A2","SLC7A3") ] <- "SLC7A*"
#sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL %in% c("IMPDH1","IMPDH2") ] <- "IMPDH*"
sko_df_genes_ <- unique(sko_df_genes[sko_df_genes$TYPE2!='Double Drug Deletion', c('Drugs','SYMBOL_',"Subtype",'n_drugs',
                                                                                   "Predicted_Essential",
                                                                                   'Drug_rank','Drugs_one_ch')])
# # Make font map for the essential genes
# bold_map <- sko_df_genes_[,]  |> 
#   group_by(SYMBOL_,Predicted_Essential)|>
#   mutate(n_drugs=n_distinct(Drugs)) |>
#   select(SYMBOL_,Predicted_Essential,n_drugs)|>
#   distinct() |>
#   dplyr::arrange(desc(n_drugs))
# bold_map$Predicted_Essential <- factor(bold_map$Predicted_Essential,c('Essential drug targets','Non-essential drug targets'))
# 
# #bold_map <- unique(bold_map[,c("SYMBOL_","n_drugs","Predicted_Essential")])
# bold_map$FACE <- ""
# bold_map$FACE[bold_map$Predicted_Essential=='Essential drug targets'] <- "bold"
# bold_map$FACE[bold_map$Predicted_Essential=='Non-essential drug targets'] <- "plain"
# bold_map$FACE <- factor(bold_map$FACE,c('bold','plain'))
n_distinct(sko_df_genes_$SYMBOL_)
unique(sko_df_genes_$SYMBOL_)
unique(sko_df_genes$SYMBOL[sko_df_genes$TYPE2 =="Single Drug Deletion"])
n_distinct(sko_df_genes$SYMBOL[sko_df_genes$TYPE2 =="Single Drug Deletion"])

#bold_map = data.frame(FACE = rep("plain",n_distinct(sko_df_genes_$SYMBOL_)))
bold_map <-sko_df_genes_# unique(sko_df_genes_[,c("SYMBOL_","Subtype","n_drugs","Predicted_Essential")])
#bold_map %>% group_by(SYMBOL_) %>% mutate(n_drugs=max(n_drugs)) -> bold_map
#bold_map <- unique(bold_map)
bold_map$FACE <-"plain"

bold_map$FACE[bold_map$Predicted_Essential=='Essential drug targets']  <- "bold"
bold_map <- bold_map[order(bold_map$FACE),]
n_essential <- n_distinct(sko_df_genes_$SYMBOL_[sko_df_genes_$Predicted_Essential=='Essential drug targets'])
bold_map$FACE[1:n_essential] <- "bold"
bold_map$FACE <- factor(bold_map$FACE,c('bold','plain'))
# 
# # Make font map for the essential genes
# sko_df_genes_bold <- unique(sko_df_genes_[,c("SYMBOL_","n_drugs","Predicted_Essential")])
# sko_df_genes_bold$Predicted_Essential <- factor(sko_df_genes_bold$Predicted_Essential,
#                                                 c('Essential drug targets','Non-essential drug targets'))
# sko_df_genes_bold <- sko_df_genes_bold[,]  |> 
#   group_by(SYMBOL_) |>
#   mutate(n_drugs=max(n_drugs)) |>
#   select(SYMBOL_,n_drugs,Predicted_Essential)|>
#   distinct()# |>
# #fct_reorder2(SYMBOL_,Predicted_Essential,desc(n_drugs))
# 
# sko_df_genes_bold$Essential <-1
# sko_df_genes_bold$Essential[sko_df_genes_bold$Predicted_Essential=='Non-essential drug targets'] <-0
# 
# fct_reorder2(sko_df_genes_bold$SYMBOL_,sko_df_genes_bold$n_drugs,sko_df_genes_bold$Essential,.desc = T)
# 
# 
# sko_df_genes_bold <- sko_df_genes_bold[fct_reorder2(,n_drugs),]
# sko_df_genes_bold$FACE <- ""
# sko_df_genes_bold$FACE[sko_df_genes_bold$Predicted_Essential=='Yes'] <- "bold"
# sko_df_genes_bold$FACE[sko_df_genes_bold$Predicted_Essential=='No'] <- "plain"
# 
sko_df_genes_$Subtype_ <-str_c("i",sko_df_genes_$Subtype)
sko_df_genes_$Subtype_ <- factor(sko_df_genes_$Subtype_,c("iGBM","iAST","iODG"))

P2_2 <- ggplot(sko_df_genes_,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
               aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank))+
  geom_tile(alpha=0.7,stat = 'identity',aes(fill=Drugs)) +
  geom_text(aes(label=Drugs_one_ch,color=Drugs,group=Drugs_one_ch, size=Predicted_Essential),position="identity" )+ #,size=2.5
  facet_grid(Predicted_Essential~Subtype_,scales = 'free')+#, space='free_y')+
  
  force_panelsizes(rows = c(1.5, 2.7),
                   cols = c(1, 1)) +

  ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
  ylab("")+xlab("Number of drugs per gene")+
  theme_classic() +
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=rep("black",33),labels=sort(x$Drugs))+
  scale_size_manual(values=c(4,2.5))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 9) ,
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12))+
        #axis.text.y= element_markdown(size = 8,face =bold_map$FACE )
        #axis.text.y = element_text(angle = 0,size = 8,face =bold_map$FACE )
        #axis.text.y = element_text(angle = 0, 
        #                           color = ifelse(sko_df_genes_$defensive_industries == "N", "red", "black"))

  guides(fill=guide_legend(title="Drugs",ncol = 1,keywidth = 0.5),
         color=guide_legend(override.aes = list(label = sort(x$Drugs_one_ch))),size="none"
         ) 
P2_2
#fill <- c("#FF0000","#2FD4DC","#FF9300","#0000FF","#00B900","#FF30F7","#6100A0","#AD7E55")
fill <- c("#2FD4DC","#00B900","#FF30F7","#FF0000","#2FD4DC","#FF0000","#FF0000","#FF0000","#2FD4DC","#FF0000",
          "#FF0000","#FF9300","#FF0000","#FF30F7","#0000FF","#FF0000","#FF0000","#00B900","#00B900","#2FD4DC",
          "#0000FF","#00B900","#00B900","#6100A0","#FF0000","#00B900","#2FD4DC","#2FD4DC","#2FD4DC","#AD7E55",
          "#2FD4DC","#2FD4DC","#2FD4DC")
unique(sko_df_genes_$Drugs)
legend_text_colors <- list()
legend_text_colors["ASC"] <- "red"
legend_text_colors["ASC2"] <- "blue"

drugs <-  unique(sko_df_genes_$Drugs)
# legend_text_colors <- data.frame(drugs,fill)
# legend_text_colors <- lapply(1:nrow(legend_text_colors), function(i) {
#   row_dict <- as.list(df[i, ])
#   names(row_dict) <- colnames(df)
#   row_dict
# })
# 
# colnames(legend_text_colors) <-
# legend_text_colors <- as.list(legend_text_colors)
# legend_text_colors[1]
# legend_text_colors <- list(
#   "A" = "purple",
#   "B" = "orange",
#   "C" = "darkblue"
# )

#P2_2 + theme(legend.text=element_markdown(color=legend_text_colors))

#P2_2 <- P2_2 + theme(legend.text=element_text(color=fill))

#g <- ggplotGrob(P2_2)
#g$grobs[[25]]
#sko_df_genes$SYMBOL_ <- sko_df_genes$SYMBOL
#sko_df_genes$SYMBOL_[str_detect(sko_df_genes$SYMBOL_,"^CA.*[0-9].*")] <- "CA genes"

x2 <- unique(sko_df_genes[sko_df_genes$TYPE2=='Double Drug Deletion',
                          c("Drugs","Drugs_one_ch")])
#sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL_ %in% c("CA*","CA**","CA13") ] <- "CA*"
sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL %in% c("SLCO1A2","SLCO1B1",'SLCO2A1','SLCO2B1') ] <- "SLCO*"

sko_df_genes_2 <- unique(sko_df_genes[sko_df_genes$TYPE2=='Double Drug Deletion', c('Drugs','SYMBOL_','Subtype','Drugs_one_ch')])

sko_df_genes_2$Subtype_ <-str_c("i",sko_df_genes_2$Subtype)
sko_df_genes_2 %>% group_by(SYMBOL_,Subtype) %>%
  mutate(
    SYMBOL_id = str_length(paste0(unique(Drugs),collapse = "; " )),
    Drug_len = str_length(Drugs),n_drugs = n_distinct(Drugs),
    Drug_rank = ntile(rank(Drugs),n_distinct(Drugs)),
  ) ->sko_df_genes_2

sko_df_genes_2$Subtype_ <- factor(sko_df_genes_2$Subtype_,c("iGBM","iAST","iODG"))
sko_df_genes_2 <- sko_df_genes_2[,!colnames(sko_df_genes_2)%in% "SYMBOL"]
sko_df_genes_2 <- unique(sko_df_genes_2)

P3_2 <- ggplot(sko_df_genes_2,
               aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank))+
  geom_tile(alpha=0.7,stat = 'identity',aes(fill=Drugs)) +
  geom_text(aes(label=Drugs_one_ch,color=Drugs,group=Drugs_one_ch),position="identity" ,size=3.7)+
  ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
  ylab("")+xlab("Number of drugs per gene")+
  facet_grid(.~Subtype_,scales = 'free')+#, space='free_y')+
  theme_classic() +
  scale_fill_manual(values=mycols[1:19])+
  scale_color_manual(values=rep("black",19),labels=sort(x2$Drugs))+
  #scale_fill_discrete(labels=c('High Program', 'Low Program'))
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 11),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=10.5),legend.title = element_text(size=12),
        axis.title = element_text(size = 12))+
  guides(fill=guide_legend(title="Combination\ndrugs",ncol = 1,keywidth = 0.5),#,override.aes = list(fill=NA)
         color=guide_legend(title="Combination\ndrugs",override.aes = list(label = sort(x2$Drugs_one_ch)))
  ) + #theme(legend.position = c(0.54,0.33))
  theme(legend.position = c(0.2,0.5))
P3_2

# 
# P3_2 <- ggplot(sko_df_genes_,
#                aes(y=reorder(fct_reorder(SYMBOL_,Drugs),n_drugs),x=fct_reorder(Drugs,SYMBOL_),
#                    group=Drugs,fill=Subtype))+
#   geom_tile(alpha=0.7,position = position_jitterdodge(reverse = F,preserve = "total")) +
#   #geom_text(aes(label=Drugs_one_ch))+
#   #facet_grid(.~Subtype,scales = 'free')+
#   ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
#   ylab("Targets")+xlab("")+#Number of drugs per gene
#   scale_fill_manual(values = c('#BB173A','#2F6790','forestgreen'))+#,'green4'))+
#   theme_classic() +
#   #scale_fill_manual(values=mycols)+
#   theme(axis.text.x = element_text(angle = 90,size = 11),
#         axis.text.y = element_text(angle = 0,size = 12) )+
#   coord_flip()
# #theme(legend.position = c(0.7,0.5),legend.background = element_blank(),legend.key = element_blank())
# 
# #guides(fill="none")
# P3_2 

#P2+ P2_2 + plot_layout(guides = "collect",ncol = 2,)
# v1 <- c('Date', 'vix1', 'vix2', 'vix3', 'doSG124', 'doSG220','vix56')
# v1 <- c("CA8","CA8A","CA7")
# v1 <- c("BRCA","BRCB","CA7")
# 
# tbl <- table(sub('\\d+$', '', v1))
# names(which.max(tbl))
# #sko_df_genes$SYMBOL_[sko_df_genes$SYMBOL=='ACAT1'] <- 'ACAT1'
# 
# P3_2 <- ggplot(sko_df_genes[sko_df_genes$TYPE2=='Double Drug Deletion', ],
#                aes(x=reorder(fct_reorder(SYMBOL_,Drugs),n_drugs),y=fct_reorder(Drugs,SYMBOL_),
#                    group=Drugs,fill=Subtype))+
#   geom_tile(alpha=1,stat = 'identity') +
#   #geom_text(aes(label=Drugs_one_ch))+
#   facet_grid(.~Subtype,scales = 'free')+
#   ggtitle("B") +#  #ggtitle("Predicted drug targets") +#
#   xlab("Targets")+ylab("")+#Number of drugs per gene
#   scale_fill_manual(values = c('#BB173A','#2F6790','forestgreen'))+#,'green4'))+
#   theme_classic() +
#   #scale_fill_manual(values=mycols)+
#   theme(axis.text.x = element_text(angle = 90,size = 11),
#         axis.text.y = element_text(angle = 0,size = 12) )+
#   guides(fill="none")  #guide_legend(title="Drugs",ncol = 1
# P3_2

# P2_final <-P2 + #xlim(0, 1300) +
#   geom_text_repel(aes(label=Evidence2),size=3.5,alpha=0.9,color='black',
#                   #force        = 4,
#                   nudge_x      = 3,hjust=0,
#                   nudge_y      = 0.2,
#                   direction    = "both",
#                   #hjust        = 0,
#                   fontface = 'bold',
#                   segment.colour = "black",
#                   box.padding = unit(0.35, "lines"),
#                   segment.size = 0.2
#   ) 

## Essential genes

## Bar plot of the number of predicted esssential genes ##
data_longer_[,] %>% #c("Model", "Data","Curation","Subtype")
  group_by(Subtype, Model, Data,Curation) %>%
  filter( Predicted_Essential==1) %>%
  dplyr::summarise(All_Genes=n_distinct(Gene_ENTREZ),
                   Common_Essential = n_distinct(Gene_ENTREZ[Essentiality_Type=='Common-Essential']))-> data_longer_genes

write_csv(data_longer_[data_longer_$Predicted_Essential==1,],'Integrated_Drug_Db/Predicted_essential_genes.csv')

data_longer_genes <- melt(data_longer_genes)

p <- ggplot(data_longer_genes, aes(x=Subtype,y=value,fill = variable)) + #,x = (..count..)/sum(..count..))
  geom_bar(position = "dodge", stat="identity")+
  geom_text(aes(label=value), position=position_dodge(width=0.9), stat="identity", vjust=-0.5)+
  ggtitle("Predicted Essential Genes in glioma subtypes") +
  #facet_grid(Model ~Data ,scales="free") +
  ylab("Number of predicted essential genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 14))+
  theme(legend.title = element_blank()) 
p
## Heat map of genes vs subtype
data_longer_$Essentiality_Type[data_longer_$Essentiality_Type =='NA_'] <- 'Normal'
data_longer_$Essentiality_Type[data_longer_$Essentiality_Type =='Normal'] <- 'Cancer type-specific'
data_longer_$Essentiality_Type[data_longer_$Essentiality_Type =='Common-Essential'] <- 'Commone essential'
data_longer_$Essentiality_Type[data_longer_$Essentiality_Type =='Commone essential'] <- 'Common essential'
data_longer_$Subtype <- factor(data_longer_$Subtype ,c('GBM','AST','ODG'))
data_longer_$Symbol <- str_replace(data_longer_$Symbol," ","")
## Merge literature Information to the essential genes
genes_literature <- na.omit(genes_literature)
genes_literature$Subtype <- 'ODG'
colnames(genes_literature)[1] <- 'Symbol'
colnames(genes_literature)[4] <- 'Evidence'

data_longer_ <- left_join(data_longer_,genes_literature[c('Symbol','Evidence')])
data_longer_$Common_Essential <- 'No'
data_longer_$Common_Essential[data_longer_$Essentiality_Type =='Common essential'] <- 'Yes'
data_longer_$Common_Essential <- factor(data_longer_$Common_Essential ,c('Yes','No'))

data_longer_$Glioma_Essential <- 'No'
data_longer_$Glioma_Essential[!is.na(data_longer_$Evidence)] <- 'Yes'
data_longer_$Glioma_Essential <- factor(data_longer_$Glioma_Essential ,c('Yes','No'))
data_longer_$Evidence2 <- str_wrap(data_longer_$Evidence , width = 40)

data_longer_$Evidence2[data_longer_$Subtype!='ODG' & data_longer_$Symbol !='SLC6A14'] <- NA
#data_longer_$Evidence2[data_longer_$Symbol %in% c('SPTLC2',"SPTLC3")] <- data_longer_$Evidence2[data_longer_$Symbol=='SPTLC1'] 
c(data_longer_$Symbol[data_longer_$Predicted_Essential==1])
data_longer_ <- left_join(data_longer_,sko_df_genes[,c('SYMBOL','n_drugs')],
                          by=c("Symbol"="SYMBOL"))
data_longer_$n_drugs[is.na(data_longer_$n_drugs)] <- 0
data_longer_$Glioma_Essential_2 <- 'No evidence'
data_longer_$Glioma_Essential_2[data_longer_$Glioma_Essential == 'Yes'] <- 'Decrease proliferation'
data_longer_$Glioma_Essential_2[data_longer_$Symbol == 'PCYT2'] <- 'Increase proliferation'

font_import()
fonts()
bold_map <- rep("plain",n_distinct(data_longer_$Symbol[data_longer_$Predicted_Essential==1]))
#n_essential <- n_distinct(sko_df_genes_$SYMBOL_[sko_df_genes_$Predicted_Essential=='Essential drug targets'])
bold_map[1:n_essential] <- "bold"

data_longer_$Subtype_ <-str_c("i",data_longer_$Subtype)
data_longer_$Subtype_ <- factor(data_longer_$Subtype_,c("iGBM","iAST","iODG"))

P2 <- ggplot(data_longer_[data_longer_$Predicted_Essential==1,], 
             aes(x=Subtype_,y=reorder(Symbol,n_drugs),color=Subtype,
                 shape=Glioma_Essential_2,
                 #alpha=Glioma_Essential,
                 group=Evidence2,
                 size = Common_Essential)) +
  geom_point()+ # position = position_dodge(width = 0.5)
  scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'),guide='none')+
  scale_size_manual(values = c(10,6))+
  scale_shape_manual(values=c("\u2193","\u2191","\u2022" ))+ # "\u23EC","\u23EB","\u2022" "\uf0e3","\uf0e4",
  scale_alpha_manual(values = c(1,0.3))+
  #scale_fill_brewer(palette="Dark2")+
  ggtitle("A") +ylab('Predicted essential genes')+ xlab('Subtype models')+#Predicted essential genes
  #facet_grid(Model ~Data ,scales="free") +
  #ylab("Predicted essential genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12))+
        #axis.text.y = element_markdown(angle = 0,size = 12, face =rev(bold_map)),
  guides(fill=guide_legend(title="Gene essentiality"),
         size=guide_legend(title="Common essential",title.position = "top",ncol = 1,
                           override.aes = list(shape="\u2022")),#"\u2195"
         #alpha=guide_legend(title="Essential in glioma\nfrom literature")
         shape=guide_legend(title="Knock-out/down\nin glioma",ncol=1
                            ,title.position = "top",override.aes = list(size=8)))+
theme(legend.position='bottom',legend.justification = 0.8)


#P2_final <- P2 + plot_spacer()+ P2_2+ plot_layout(ncol = 3,tag_level = "new",#guides = "collect",                                                  widths = c(3.5, -0.75 ,4.8))#
#P2_final
P2_final <- cowplot::plot_grid(P2, P2_2, ncol = 2,rel_widths = c(3.3,4.8))
P2_final
ggsave(P2_final,filename="Figure/Figure_2_Essential_Genes.png", units="in", width=11, height=10.5,dpi=300)
# quartz(type = 'pdf', file = 'Figure/Figure_2_Essential_Genes_2.pdf', width=11, height=10.5)
# P2_final
# dev.off()

# ggsave(P2_final,filename="Figure/Figure_2_Essential_Genes.pdf", units="in", width=11, height=11.5,dpi=300)
# library(Cairo)
# 
# quartz(type = 'pdf', file = 'Figure/Figure_2_Essential_Genes.pdf', width=11, height=11.5)
# P2_final
# dev.off()


### Figure 3 
# Predicted combination are synergistic with drug targets
dko_df  = read_csv('Integrated_Drug_Db/DKO_Combination_result.csv')
colnames(dko_df) <- c('Drugs','Subtype','grRatio')
dko_df$Drug1 <- str_split(dko_df$Drugs,' ;',simplify = T)[,1]
dko_df$Drug2 <- str_split(dko_df$Drugs,' ;',simplify = T)[,2]
unique(dko_df$Drug1)
p2 <- ggplot(dko_df,aes(x=Drug1,y=Drug2,fill=Subtype))+
  geom_tile()+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(. ~Subtype ,scales="free")+
  ggtitle("Drug combination prediction") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))
p2

## Compare the grRion of the drug combinations individually
dko_df_sko  = read_csv('Integrated_Drug_Db/DKO_Combination_individual_grRatio_result.csv')
colnames(dko_df_sko) <- c('Drugs','Subtype','grRatio_single')

# Add the dko without cutoff on grRatio
n_distinct(dko_df_sko$Drugs)
dko_df_all  = read_csv('Integrated_Drug_Db/DKO_Combination_result_all_predictions.csv')
colnames(dko_df_all) <- c('Drugs','Subtype','grRatio','DelRxns')
dko_df_all <- dko_df_all[,colnames(dko_df_all)!='DelRxns']
dko_df_all <- unique(dko_df_all[,colnames(dko_df_all) != "DelRxns"])

dko_df_all$Drug1 <- str_split(dko_df_all$Drugs,' ;',simplify = T)[,1]
dko_df_all$Drug2 <- str_split(dko_df_all$Drugs,' ;',simplify = T)[,2]
dko_df_all <- dko_df_all[dko_df_all$Drug1 %in% dko_df_sko$Drugs,]
dko_df_all <- dko_df_all[dko_df_all$Drug2 %in% dko_df_sko$Drugs,]
dko_df_all <- dko_df_all[dko_df_all$Drug1 != dko_df_all$Drug2,]
## Make sure that combinations with no measure grRatio due to missing targets for faster computing
## in the model are added
drugs_1 = unique(dko_df$Drug1)
drugs_2 = unique(dko_df$Drug2)
x <- data.frame(Drug1=paste0(drugs_1,collapse = "; "),
                Drug2=paste0(drugs_2,collapse = "; "),
                Subtype=paste0(c("AST","GBM","ODG"),collapse = "; "))
x %>%separate_rows(Drug1,sep = "; ") %>%separate_rows(Subtype,sep = "; ")%>%separate_rows(Drug2,sep = "; ")-> x
dko_df_all <- left_join(x,dko_df_all)
dko_df_all$Drugs[is.na(dko_df_all$grRatio)] <- str_c(dko_df_all$Drug1[is.na(dko_df_all$grRatio)] ,dko_df_all$Drug2[is.na(dko_df_all$grRatio)] ,sep = " ;")
dko_df_all$grRatio[is.na(dko_df_all$grRatio)] <- 1

#drugs_all = unique(dko_df_sko$Drugs)3
#x <- data.frame(Drugs=paste0(drugs_all,collapse = "; "),
#                 Subtype=paste0(c("AST","GBM","ODG"),collapse = "; "))
#x%>%separate_rows(Drugs,sep = "; ") %>%separate_rows(Subtype,sep = "; ")-> x
# 
# dko_df_sko <- left_join(x,dko_df_sko)
# dko_df_sko$Drugs[is.na(dko_df_all$grRatio)] <- dko_df_sko$Drug1[is.na(dko_df_sko$grRatio)] ,dko_df_all$Drug2[is.na(dko_df_all$grRatio)] ,collapse = " ;")
# dko_df_all$grRatio[is.na(dko_df_all$grRatio)] <- 1
# 

dko_df_all <- data.frame(dko_df_all)
nrow(dko_df_all)

dko_df <- data.frame(dko_df)

dko_df_sko <- data.frame(dko_df_sko)
dko_df_ <-left_join(dko_df_all,dko_df_sko,by=c('Drug1'='Drugs','Subtype'='Subtype'))

dko_df_ <-left_join(dko_df_,dko_df_sko,by=c('Drug2'='Drugs','Subtype'='Subtype'))
colnames(dko_df_)

dko_df_ <- dko_df_[,c("Drugs" , "Subtype"   ,  "grRatio" , "Drug1" ,"Drug2","grRatio_single.x" ,"grRatio_single.y") ]
dko_df_ %>%  pivot_longer(cols =c("grRatio","grRatio_single.x",'grRatio_single.y') ,
                          names_to = "Deletion_type",values_to = "grRatio") ->dko_df_longer
dko_df_longer$Deletion_type[dko_df_longer$Deletion_type=='grRatio'] <- 'Combination'
dko_df_longer$Deletion_type[dko_df_longer$Deletion_type%in%c("grRatio_single.x")] <- 'Drug1'
dko_df_longer$Deletion_type[dko_df_longer$Deletion_type%in%c("grRatio_single.y")] <- 'Drug2'
dko_df_longer$Subtype <- factor(dko_df_longer$Subtype,c('GBM','AST','ODG'))

dko_df_longer <- dko_df_longer[dko_df_longer$Drugs %in% dko_df$Drugs ,]
dko_df_longer$Drug1 <- firstup(dko_df_longer$Drug1)
dko_df_longer$Drug2 <- firstup(dko_df_longer$Drug2)
dko_df_longer$Combination_name <- str_c(dko_df_longer$Drug1,str_to_lower(dko_df_longer$Drug2),sep  = "/")
dko_df_longer$growth_red <- 1 - dko_df_longer$grRatio
#dko_df_longer$EA <- dko_df_longer$growth_red[dko_df_longer$Deletion_type == 'Drug1'] 
            
dko_df_longer %>% 
 # mutate(  ) %>%
  group_by(Combination_name, Subtype) %>% mutate(
    EA = unique(growth_red[Deletion_type == 'Drug1']),
    EB = unique(growth_red[Deletion_type == 'Drug2']),
    EAB = unique(growth_red[Deletion_type == 'Combination']),
  #Combination Subthresholding
  Synergism_score = sum(grRatio[Deletion_type == 'Combination'])- sum(grRatio[Deletion_type == 'Drug1'] , grRatio[Deletion_type == 'Drug2']),
  CI_Subthresholding = (EA + EB)/EAB,
  #Bliss Independence
  Sum_A_B = EA+EB/EAB,
  CI_Bliss  = (EA+EB-(EA*EB))/EAB,
  ) -> dko_df_longer
dko_df_longer_result <- dko_df_longer[dko_df_longer$Deletion_type == 'Combination',]
df2 <- unique(dko_df_longer[order(dko_df_longer$Synergism_score),c('Combination_name','Drug1','Drug2')])
myColors <- data.frame(myColors=c("#21de04","#2bdba6","#ffcc33", "#ed9b0e" ,"#cf68ed"),
                       Gene_type=c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
                                   "Approved anti-melanoma target","NO-related genes"))

dko_df_longer$Combination_name
dko_df_longer$Combination_name_ <- str_replace(dko_df_longer$Combination_name," \\+ ","\\/")
dko_df_longer$Subtype_ <-str_c("i",dko_df_longer$Subtype)
dko_df_longer$Subtype_ <- factor(dko_df_longer$Subtype_,c("iGBM","iAST","iODG"))

#dko_df_longer$Synergism_score[dko_df_longer$Combination_name =="Fluorouracil + Celecoxib"] <- -3.97610
# Rank the combinations based on the selection done
comp_selection <- read_excel("Supplementary File 2.xlsx",sheet = 4)
dko_df_longer$Combination_name <- str_c(dko_df_longer$Drug1,"-",str_to_lower(dko_df_longer$Drug2))
dko_df_longer <- left_join(dko_df_longer,comp_selection[,c("Combinations","Rank")], by=c("Combination_name"="Combinations"))
dko_df_longer$Deletion_type <- factor(dko_df_longer$Deletion_type,c("Drug1","Drug2","Combination" ))
P3_1 <- ggplot(dko_df_longer,
       aes(x=(1- grRatio) *100,color= forcats::fct_rev(Deletion_type)
           ,fill= forcats::fct_rev(Deletion_type),
           y=reorder(Combination_name_,-Rank)))+  #Synergism_score # [dko_df_longer$Deletion_type!='Combination',]
  geom_bar(stat = "identity", position=position_dodge(),size=1,width = 0.8)+
  #geom_segment(aes(x=-0.05,xend=1- grRatio,y=paste0(Combination_name,Deletion_type),
  #                 yend=paste0(Combination_name,Deletion_type)),size=3.5,
  #              position=position_dodge())+
  #geom_point(aes(shape= Drug1), position=position_dodge(),size=1,width = 0.8)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(. ~Subtype_ ,scales="free")+
  xlab("Growth reduction (%)") +
  #scale_x_continuous(trans = 'log10')+
  ylab("Predicted combinations") +
  #ggtitle("The biomass growth rate in the drug combination\nprediction compared to single drug") +
  ggtitle("A")+
  #scale_y_discrete(labels= df2$Drug2)+
  theme_classic()+
  #scale_x_continuous(trans = 'log10')+
  theme(axis.text.x = element_text(angle = 0,size = 10),#legend.text = element_text(angle = 0,size = 10),
        axis.text.y = element_text(angle = 0,size = 12),axis.title.x =element_text(angle = 0,size = 12) ,
        strip.text = element_text(size = 12),
        legend.text=element_text(size=12),legend.title = element_text(size=12),
        axis.title = element_text(size = 12))+
  guides(color=guide_legend(title = 'Drug prediction',ncol=3,title.position = "left",reverse = TRUE),
         fill=guide_legend(title = 'Drug prediction',ncol=3,title.position = "left",reverse = TRUE))+
  theme(plot.margin=unit(c(0.2,0.21,1,0.2),"cm")) +
  theme(legend.position=c(-0.01,-0.08),legend.box = "horizontal") #legend.position=c(0.05,-0.15)
  
  #expand_limits(x = -0.2, y = 0)
  #scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1)) 
P3_1
empty <- theme(
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
# P3_Final <- P3_1 + plot_spacer()+ P3_2+ plot_layout(ncol = 3,tag_level = "new",#guides = "collect",
#                                                widths = c(4, -0.1 ,3))#
# P3_Final
P3_Final <- cowplot::plot_grid(P3_1, P3_2+ theme(legend.position = c(0.21,0.5)), ncol = 2,rel_widths = c(4 ,4.4))#,-.013 #+ theme(legend.position = c(0.2,0.5))
#P3_Final
ggsave(P3_Final,filename="Figure/Figure_3_Synergistic_Combinations.png", units="in", width=11.5, height=8.5,dpi=300,)
#ggsave(P3_Final,filename="Figure/Figure_3_Synergistic_Combinations_2.pdf", units="in", width=11.5, height=8.5,dpi=300,)


### Figure 4
# Scatter plot of drug potency and CSF bioavailability
### Point plot of the preclinical data with DOSAGE  in x axis and Color to EFFECT and shape to the top frequent cell lines
vitro <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S14_inVitro_data")
vitro$Drug <- firstup(vitro$Drug)

vitro %>% 
  filter(!is.na(Dosage_in_µM)) ->vitro
colnames(vitro)
vitro <- vitro[,c('Drug','Cell_lines','Dosage_in_µM','Effect','Prediction','Effect with TMZ')]
vitro$Cell_lines[vitro$Cell_lines=="patient-derived GSC-enriched 3D spheroid cultures"] <- 'GSC spheroid'
vitro$Dosage_in_µM <- str_replace(vitro$Dosage_in_µM,",",".")

vitro$Dosage_in_µM <- as.numeric(vitro$Dosage_in_µM)

### Classify cell lines based on the glioma subtype
depmap_cellinfo <- read.csv("CRISPR/DepMap_22Q1/sample_info_braincancer.csv")
vitro <- left_join(vitro,depmap_cellinfo,by=c('Cell_lines'='cell_line_name'))

vitro$Subtype[str_detect(vitro$Cell_lines,'GSC')] <- 'GBM'


unique(vitro$lineage_subtype)

nrow(vitro[is.na(vitro$Subtype),])
nrow(vitro[!is.na(vitro$Subtype),])
colnames(vitro)

colnames(vitro)[6] <- "Effect_with_TMZ"
## Select combination classes most frequent
x <- as.data.frame(table(vitro$Cell_lines))
x <- x[x$Var1 !='-',]

x_list <- x[x$Freq<3,'Var1']

vitro$Cell_Lines <- vitro$Cell_lines
vitro$Cell_Lines[vitro$Cell_lines %in% x_list] <- 'Others'
vitro$Cell_Lines[vitro$Cell_lines == '-'] <- 'Others'
vitro$Reduction_in_Viability <- vitro$Effect
vitro_viability <- vitro[grepl("[.]",vitro$Effect),] ## Keep a placeholder to merge viability reduction with primary PRISM and Nam et al

vitro$Effect_Original <- vitro$Effect

vitro$Effect[str_detect(vitro$Reduction_in_Viability,"0.")] <- 'Reduction in Viability'


vitro$Reduction_in_Viability[str_detect(vitro$Reduction_in_Viability,'50')] <- 0.5
vitro$Effect_with_TMZ[is.na(vitro$Effect_with_TMZ)] <- 'Not tested'

vitro$Reduction_in_Viability[!str_detect(vitro$Reduction_in_Viability,'0')] <- 1
vitro$Reduction_in_Viability <- as.numeric(vitro$Reduction_in_Viability)

colnames(vitro)

vitro$Effect[vitro$Effect %in% c('Apoptosis', 'Autophagy')] <- 'Apoptosis/Autophagy'
unique(vitro$Effect)
vitro$Effect <-factor(vitro$Effect ,c("Apoptosis/Autophagy","Reduction in Viability","MIC","Increased radiosensitivity","Reducing acidosis and hypoxia",'No effect',"Proliferation" ,'Not tested'))

unique(vitro$Effect)
list_1 <- c("Reduction in Viability","MIC" ,"Apoptosis/Autophagy")

length(unique(vitro$Drug))
# Scatterplot for the number of cell lines VS median dosage, colored by the effect of the largest dosage

vitro %>% group_by(Drug) %>%
  mutate(median_dosage = median(Dosage_in_µM), #[Effect %in% list_1]
         #Effect_of_top_dosage=if_else(length(unique(Effect)) >1,Effect[max(Dosage_in_µM)],Effect[1]),
         Effect_of_top_dosage=Effect[which.max(Dosage_in_µM)],
         n_celllines=n_distinct(Cell_lines))->vitro_summary  #Effect %in% list_1]
vitro_summary$n_celllines [vitro_summary$Cell_lines =='-']<-0

vitro_summary <- distinct(vitro_summary[,c('Drug','median_dosage','n_celllines','Effect_of_top_dosage')])

#########################################
vitro %>% group_by(Drug) %>%
  mutate(median_dosage = median(Dosage_in_µM[Effect %in% list_1]), #[Effect %in% list_1]
         n_celllines=n_distinct(Cell_lines[Effect %in% list_1]))->vitro  #Effect %in% list_1]

### Single drugs cell lines
vitro <- as.data.frame(vitro)
#vitro[vitro$Cell_lines =='-','n_celllines'] <-0
max(vitro$median_dosage)
vitro$median_dosage[vitro$Effect =='Not tested'] <- 20000
vitro$median_dosage[vitro$Effect =='Proliferation'] <-   15000
vitro$median_dosage[vitro$Drug %in% c('levodopa')] <- 13000
vitro$median_dosage[vitro$Drug %in% c("gabapentin","inositol","zalcitabine")] <- 12000
#vitro[vitro$Drug %in% c("didox"),'median_dosage'] <- 10000
#vitro[vitro$Drug %in% c("levodoa",'hydroxurea'),'median_dosage'] <- 10000
vitro$median_dosage[vitro$Drug %in% c("decitabine","acyclovir")] <- 8000
vitro$median_dosage[vitro$Drug %in% c("didox","hydroxyurea")] <- 7500

vitro$median_dosage[vitro$Drug %in% c("gemcitabine","melphalan",'arsenic-trioxide')] <- 0.001
vitro$median_dosage[vitro$Drug %in% c('ribavirin')] <- 0.002
#vitro$median_dosage[vitro$Drug %in% c('levodopa')] <- 0.003

vitro$median_dosage[vitro$Drug %in% c('fludarabine-phosphate')] <- 7000

vitro$n_celllines[vitro$Effect =='No effect'] <- 0
vitro$median_dosage[is.na(vitro$median_dosage)] <- 0
vitro$median_dosage <- as.numeric(vitro$median_dosage)
P_sko_lit <- ggplot(vitro[vitro$Prediction=='Single',],
                    aes(x=Dosage_in_µM,color=Effect,y=reorder(Drug,desc(median_dosage))))+#fct_reorder2(Drug,desc(median_dosage),desc(n_celllines))))+
  geom_point(aes(shape=Cell_Lines,size=Reduction_in_Viability),position = position_dodge(width = 0.5))+
  ylab("Drugs") +
  xlab("Dosage in µM") +
  #facet_wrap(.~Prediction,scales = 'free')+
  scale_color_manual(values = c("chartreuse4", "green3", "olivedrab2",'red','darkred','grey'))+
  scale_shape_manual(values = seq(1,n_distinct(vitro[vitro$Prediction=='Single',"Cell_Lines"]),1)) +
  scale_x_continuous(trans='log10',n.breaks = 10)+
  geom_vline(xintercept=50, linetype="dashed",  color = "red", size=1)+
  geom_vline(xintercept=10, linetype="dashed",  color = "blue", size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 13),
        axis.title.x =element_text(angle = 0,size = 12) ,
        strip.text = element_text(size = 12),
        legend.text=element_text(size=12),legend.title = element_text(size=12),
        axis.title = element_text(size = 12))+
  ggtitle("Efficacy of predicted single drugs in in vitro cell viability assay from literature")
ggsave(P_sko_lit,filename="Figure/Figure_x_Viability_singleDrugs_Literature.png", units="in", width=9, height=8,dpi=300,)

# vitro_nano <- vitro[vitro$Prediction=='Single' & vitro$Dosage_in_µM<=0.2& vitro$Dosage_in_µM >0,]
# ggplot(vitro_nano,aes(x=Dosage_in_µM,color=Effect,y=reorder(Drug,desc(median_dosage))))+#fct_reorder2(Drug,desc(median_dosage),desc(n_celllines))))+
#   geom_point(aes(shape=Cell_Lines,size=Reduction_in_Viability),position = position_dodge(width = 0.5))+
#   ylab("Drugs") +
#   xlab("Dosage in µM") +
#   #facet_wrap(.~Prediction,scales = 'free')+
#   scale_color_manual(values = c("chartreuse4",'red','darkred'))+
#   scale_shape_manual(values = seq(1,22,1)) +
#   scale_size_continuous(range=c(3,5))+
#   #scale_x_continuous(trans='log10',n.breaks = 10)+
#   #geom_vline(xintercept=50, linetype="dashed",  color = "red", size=1)+
#   #geom_vline(xintercept=10, linetype="dashed",  color = "blue", size=1)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 13))+
#   ggtitle("Efficacy of predicted single drugs in in vitro cell viability assay from literature")
# 
# unique(vitro[vitro$Prediction=='Single','Cell_Lines'])
# #vitro_s <- distinct(vitro[,c('Drug','median_dosage')])
# #vitro_s <- vitro_s[order(vitro_s$median_dosage),]
# #write_clip(vitro_s)

## combination cell lines
x <- vitro[vitro$Prediction=='Combination',]
x$median_dosage[x$Effect =='Not tested'] <- 20000
x$median_dosage[x$Drug%in% c("cannabidiol" , "fluorouracil" )] <- 0.001
x$median_dosage[x$Drug%in% c("eflornithine" )] <- 0.002
x$median_dosage[x$Drug%in% c("resveratrol" )] <- 0.003
x$median_dosage[x$Drug%in% c("brinzolamide")] <-19000
x$median_dosage[x$Drug%in% c("zidovudine","acetazolamide")] <-18000


max(x$median_dosage)
unique(x$Drug)
P_dko_lit <- ggplot(x,aes(x=Dosage_in_µM,color=Effect,y=reorder(Drug,desc(median_dosage))))+#fct_reorder2(Drug,desc(median_dosage),desc(n_celllines))))+
  geom_point(aes(shape=Cell_Lines,size=Reduction_in_Viability),position = position_dodge(width = 0.6))+
  ylab("Drugs") +
  xlab("Dosage in µM") +
  #facet_wrap(.~Prediction,scales = 'free')+
  scale_color_manual(values = c("chartreuse4", "green3", "olivedrab2",'deepskyblue4','deepskyblue1','red','grey'))+
  #theme_classic() +
  scale_shape_manual(values = seq(1,23,1)) +
  scale_x_continuous(trans='log10',n.breaks = 10)+
  geom_vline(xintercept=50, linetype="dashed",  color = "red", size=1)+
  geom_vline(xintercept=10, linetype="dashed",  color = "blue", size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 13),
        axis.text.y = element_text(angle = 0,size = 12))+
  ggtitle("Efficacy of predicted drug combination in in vitro cell viability assay from literature") 
ggsave(P_dko_lit,filename="Figure/Figure_x_Viability_Combinations_Literature.png", units="in", width=9, height=8,dpi=300,)

unique(vitro[vitro$Prediction=='Combination','Cell_Lines'])
unique(vitro[vitro$Prediction=='Combination','Effect_with_TMZ'])

###  Read the merged IC50 databases
Merged_ic50 <- read_csv('./Merged_IC50_databases.csv')
Merged_ic50$Drugs <- firstup(Merged_ic50$Drugs)
Merged_ic50_drugs <- unique(Merged_ic50[,c('Drugs',"Database")])
#Merged_ic50$Drugs[Merged_ic50$Drugs =='Rifampin'] <- "Rifamycin"
Merged_ic50$Drugs[Merged_ic50$Drugs =='5-fluorouracil'] <- "Fluorouracil"
Merged_ic50$Drugs[Merged_ic50$Drugs =='Fludarabine'] <- "Fludarabine-phosphate"


bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
drugs_all <- bbb_data[,c('Drugs','Prediction')]
drugs_all$Drugs <- firstup(drugs_all$Drugs)

#Merged_ic50$Drugs <- str_to_lower(Merged_ic50$Drugs )
n_distinct(Merged_ic50$Drugs)
n_distinct(Merged_ic50$cell_line_name)

Merged_ic50$Drugs <- str_replace(Merged_ic50$Drugs,'-',' ')
table(unique(Merged_ic50[,c('Database','Drugs')])[,c('Database')]) 

drugs_ic50 <- left_join(drugs_all,Merged_ic50) # select our drugs in the merged database of IC50
n_distinct(drugs_ic50$Drugs[!is.na(drugs_ic50$IC50)])
nrow(drugs_ic50[!is.na(drugs_ic50$IC50),])

unique(drugs_ic50$Drugs[!is.na(drugs_ic50$IC50)])

# select brain cancer cell lines 
depmap_cellinfo = read.csv("./CRISPR/DepMap_22Q1/sample_info_braincancer.csv")

#depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID','cell_line_name',"Subtype",'primary_disease','CCLE_Name')]
depmap_cellinfo <- depmap_cellinfo[str_detect(depmap_cellinfo$primary_disease,'Brain Cancer'),]

depmap_cellinfo$cell_line_name <- toupper(depmap_cellinfo$cell_line_name)

drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'-','')
#drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'\\ ','')
#drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'\\.','')
drugs_ic50$cell_line_name <- toupper(drugs_ic50$cell_line_name)

#drugs_ic50 <- drugs_ic50[drugs_ic50$cell_line_name %in% depmap_cellinfo$cell_line_name,] 
drugs_ic50 <- left_join(drugs_ic50,depmap_cellinfo[c("cell_line_name",'Subtype')])
drugs_ic50 <- drugs_ic50[!is.na(drugs_ic50$Subtype),]

table(unique(drugs_ic50[,c('Database','Drugs',"Prediction")])[,c('Database',"Prediction")])
table(unique(drugs_ic50[,c('Database','Subtype',"cell_line_name")])[,c('Subtype',"Database")])

### Add literature IC50 data
vitro$Database <- 'Literature'
vitro_ic50 <- vitro[vitro$Effect_Original=='IC50',c('Drug','Prediction','Cell_lines','Dosage_in_µM','Database','Subtype')]
colnames(vitro_ic50) <- c("Drugs","Prediction" ,"cell_line_name", "IC50"  , "Database" , "Subtype")
drugs_ic50 <- rbind(drugs_ic50,vitro_ic50)

n_distinct(drugs_ic50$Drugs)
drugs_ic50 <- distinct(drugs_ic50)

drugs_ic50 %>% group_by(Drugs,Prediction) %>% 
  mutate(Median_IC50 = median(IC50[Database!="CTRPv2"],na.rm = T),
         Std_IC50 = sd(IC50[Database!="CTRPv2"],na.rm = T),
         N_cells = n_distinct(cell_line_name[Database!="CTRPv2"])) ->drugs_ic50
drugs_ic50 %>% group_by(Drugs,Prediction,Database) %>% 
  mutate(Median_IC50_per_database = median(IC50[Database!="CTRPv2"],na.rm = T)) ->drugs_ic50

drugs_ic50$Prediction <- factor(drugs_ic50$Prediction,
                                    c("Single",'Combination',"Anti-brain cancer"))
drugs_ic50$Drugs <- firstup(drugs_ic50$Drugs)
drugs_ic50$Database <- factor(drugs_ic50$Database,c("Literature","PRISM Secondary", "GDSC1000" , "GDSC2000" , "gCSI" ,"CTRPv2"  ))

P_S_IC50_1 <- ggplot(drugs_ic50[drugs_ic50$Database!= c("CTRPv2"),],#[drugs_ic50$Database %in% c("Literature","PRISM Secondary"),], 
                   aes(x=IC50, y=reorder_within(Drugs,desc(Median_IC50_per_database),Database), 
                              fill=Prediction,
                              #shape=Subtype,
                              color=Prediction)) + 
  geom_violin(alpha=0.4)+
  geom_point(alpha=0.5,size=3)+
  geom_point(aes(x=Median_IC50_per_database),color='black')+
  scale_x_continuous(trans='log10',n.breaks = 10)+
  scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  ylab("Drugs") +
  #facet_wrap(.~Database,scales = "free")+
  #force_panelsizes(rows = c(9,1.3),
  #                 cols = c(3),respect = T) +
  scale_y_reordered()+
  xlab("IC50 in µM") +theme_classic() +
  theme(axis.text.x = element_text(angle = 90,size = 10),
        axis.text.y = element_text(angle = 0,size = 12))+
  geom_vline(xintercept=1, linetype="dashed", color = "blue", size=0.5)+
  geom_vline(xintercept=10, linetype="dashed", color = "red", size=0.5)+
  ggtitle("IC50 measures in the brain cancer cell lines across four cell viability databases and literature") +
  guides(fill=guide_legend(title='Drug Type',ncol = 3),color=guide_legend(title='Drug Type',ncol = 3))+
  theme(legend.position = "bottom",legend.background = element_blank(),legend.key = element_blank())
design <- matrix(c(1,1,1,2,2,2,3,4,5), 3, 3)

P_S_IC50_1 <- P_S_IC50_1 + facet_manual(.~Database,scales = "free", design = design)

P_S_IC50_1
# P_S_IC50_2 <- ggplot(drugs_ic50[drugs_ic50$Database %in% c("GDSC1000" , "GDSC2000" , "gCSI"),], 
#                      aes(x=IC50, y=reorder_within(Drugs,desc(Median_IC50_per_database),Database), 
#                          fill=Prediction,
#                          #shape=Subtype,
#                          color=Prediction)) + 
#   geom_violin(alpha=0.3)+
#   geom_point(alpha=0.5)+
#   geom_point(aes(x=Median_IC50_per_database),color='black')+
#   scale_x_continuous(trans='log10',n.breaks = 10)+
#   scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
#   scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
#   ylab("") +
#   facet_wrap(.~Database,scales = "free",ncol = 1)+
#   scale_y_reordered()+
#   xlab("") +theme_classic() +
#   theme(axis.text.x = element_text(angle = 90,size = 10),
#         axis.text.y = element_text(angle = 0,size = 12))+
#   geom_vline(xintercept=1, linetype="dashed", color = "blue", size=0.5)+
#   geom_vline(xintercept=10, linetype="dashed", color = "red", size=0.5)+
#   ggtitle("") +
#   guides(fill=guide_legend(title='Drug Type',ncol = 1),color=guide_legend(title='Drug Type',ncol = 1))+
#   theme(legend.position = "top",legend.background = element_blank(),legend.key = element_blank())
# 
# P_S_IC50 <- cowplot::plot_grid(P_S_IC50_1, P_S_IC50_2, ncol = 2,rel_widths = c(6,3))
# P_S_IC50
#ggsave('Figure/FigureSxx_IC50.pdf', P_S_IC50,units = 'in',width = 10,height = 9,dpi = 300)
ggsave('Figure/FigureSxx_IC50.png', P_S_IC50_1,units = 'in',width = 11,height = 9,dpi = 300)

drugs_ic50 %>% group_by(Prediction,Drugs) %>% 
  summarise(Median_IC50 = median(IC50[drugs_ic50$Database!="CTRPv2"],na.rm = T),
            Std_IC50 = sd(IC50[drugs_ic50$Database!="CTRPv2"],na.rm = T),
            N_cells = n_distinct(cell_line_name[drugs_ic50$Database!="CTRPv2"])) ->drugs_ic50_summ

write_csv(drugs_ic50,'./Figure/IC50_table.csv')
write_csv(drugs_ic50_summ,'./Figure/IC50_table_summary.csv')



#write.xlsx(drugs_ic50,'./Figure/IC50_table.xlsx')
#write.xlsx(drugs_ic50_summ,'./Figure/IC50_table_summary.xlsx')


#### BBB data collected manually
bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
bbb_data$Drugs <- firstup(bbb_data$Drugs)
colnames(bbb_data)
colnames(bbb_data)[1] <- "Suggestion"
colnames(bbb_data)[4] <- "BBB_permeability"
colnames(bbb_data)[5] <- "logBB"

# vitro_data <- vitro[!vitro$Effect %in% c('Proliferation','No effect',"Increased radiosensitivity",
#                                          "Reducing acidosis and hypoxia"),]
# colnames(vitro_data)[1] <- 'Drugs'


vitro_ic50 <- drugs_ic50[,c('Drugs','cell_line_name',"Subtype" , "Prediction",'IC50','Database')] #drugs_ic50$Database!='Literature'
vitro_ic50 <- vitro_ic50[vitro_ic50$Database!= "CTRPv2",]

colnames(vitro_ic50) <- c('Drugs','Cell_lines',"Subtype" , "Prediction",'Dosage_in_µM','Database')
#vitro_ic50$Effect <- "Reduction in Viability"
#vitro_ic50$Reduction_in_Viability <- 0.5
vitro_all <- vitro_ic50# rbind(vitro_data[,colnames(vitro_ic50)],vitro_ic50)

# Drug not tested in IC50 measures
untested_drugs <- setdiff(bbb_data$Drugs,vitro_all$Drugs)#unique(vitro_all$Drugs[vitro_all$Effect=="Not tested"])

# Drug tested in IC50 databases
tested_drugs <- unique(vitro_all$Drugs) #[vitro_all$Database!="Literature"]
# 
# #Drugs untested in literature to be removed before calculating the median dosage
# removed_drugs <- intersect(untested_drugs,tested_drugs)
# removed_drugs
# idx <- rownames(vitro_all[vitro_all$Effect =="Not tested" & vitro_all$Drugs %in%removed_drugs, ])
# vitro_all <- vitro_all[!rownames(vitro_all) %in% idx,]

# take median conc for the drugs
vitro_all %>% group_by(Drugs) %>% 
  summarize(Median_Dosage =median(Dosage_in_µM,na.rm = T)) ->vitro_all_summ

bbb_data <- left_join(bbb_data[,1:5],vitro_all_summ)
bbb_data <- bbb_data[!is.na(bbb_data$Suggestion),]
bbb_data$Test_inVitro <- 'Tested'
bbb_data$Test_inVitro[bbb_data$Drugs %in% untested_drugs] <- 'Not tested' #Median_Dosage==0
bbb_data$Median_Dosage[bbb_data$Median_Dosage==0] <- 0

#

bbb_data$LogBB_Info <- 'Available'
bbb_data$LogBB_Info[bbb_data$logBB=="No info"] <- 'Not available'
bbb_data$logBB[bbb_data$logBB=="No info"] <- 0
bbb_data$logBB <- as.numeric(bbb_data$logBB)

#bbb_data <- bbb_data[bbb_data$Suggestion!='Excluded',]

bbb_data$Prediction <- factor(bbb_data$Prediction,c('Single','Combination', 'Anti-brain cancer'))

bbb_data$BBB_permeability <- factor(bbb_data$BBB_permeability,c("Permeable","Poor CSF bioavailability","Non-permeable"))#,'Approved Anti\nbrain cancer'

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
bbb_data$Neg_Log_Median_Dosage <- log10(bbb_data$Median_Dosage)
bbb_data$logBB[bbb_data$LogBB_Info =="Not available"] <- 0.5 # Drugs with no LogBB data

bbb_data$Median_Dosage[bbb_data$Test_inVitro =="Not tested"] <- 0.002 # Drugs with no IC50 data
x <- bbb_data[! (bbb_data$LogBB_Info =="Not available" & bbb_data$Test_inVitro =="Not tested") ,]

P4 <- ggplot(x,aes(x=logBB, y=Median_Dosage, color =Prediction,
                          shape= BBB_permeability,label=Drugs,size=Prediction,
                          #shape=LogBB_Info,
                          alpha=LogBB_Info))+
  geom_smooth(data =x[x$Prediction=='Anti-brain cancer',] ,
              method=lm , color="red",size=1.5,
              fill="#E9FFFF", se=TRUE,na.rm = T,aes(group=Prediction,x=logBB, y=Median_Dosage),
              fullrange=T, xseq = seq(-4.4,0.2, length=80)) +
  #coord_cartesian(xlim=c(-4.5,0), ylim=c(10000,0.01)) +
  
  #geom_smooth(data =x[x$Prediction=='Anti-brain cancer',] ,method=lm, level=0.90)+

  geom_point(size=4)+#,position = position_dodge(width = 0.4)
  
  geom_text_repel(nudge_x = 0.2,max.overlaps =53)+#nudge_x = 0.3,nudge_y = 0.25
  #ggtitle("Relationship between CSF bioavailability and median cytotoxic acivity") +
  ggtitle("") +
  xlab("CSF bioavailability: Log2( CSF conc. / Plasma conc.) ") +
  ylab("Potency: Median IC50 (µM)") +
  scale_size_manual(values=c(4,4,5))+
  geom_vline(xintercept=0.2, linetype="dashed",  color = "red", size=0.4)+
  geom_hline(yintercept=0.0075, linetype="dashed",  color = "red", size=0.4)+
  
  #geom_abline(intercept = 1.4, slope = -0.6,color='blue',linetype = "dashed") +
  scale_alpha_discrete(range = c(1,0.5))+
  scale_color_manual(values = c("forestgreen",'#ffa100','steelblue'))+#,'red'
  theme_classic() +
  scale_y_continuous(trans=reverselog_trans(base=10),expand = expansion(add = c(0.2, 0.5)),
                     labels=trans_format("identity", function(x) x),
                     breaks=c(NULL,0,0.01,0.1,1,10,100,1000,10000))+#©,
  scale_x_continuous(expand =  expansion(add = c(0.1, 0.5)),breaks=c(NULL,0,-1,-2,-3,-4))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 13) ,
        axis.title.x =element_text(angle = 0,size = 12) ,
        strip.text = element_text(size = 12),
        legend.text=element_text(size=12),legend.title = element_text(size=12),
        axis.title = element_text(size = 12))+
  annotate("text", x=0.75, y=0.04, label= "Drugs without\nLogBB data",size=5)+
  annotate("text", x=-2, y=0.0002, label= "Drugs without potency data",size=5)+
  guides(size="none",alpha=guide_legend(title="LogBB data"),
         shape=guide_legend(title="BBB permeability"))#+ ylim( 0,10000)
P4
ggsave(P4,filename="Figure/Figure_4_CSF_bioavailability_Drugs.png", units="in", width=11, height=8, dpi=300)


#### visualize the grRatio of the control model
sko_results_ctrl <- read_csv('./Integrated_Drug_Db/SKO_result_CTRL_model.csv')
sko_results_ctrl <- sko_results_ctrl[sko_results_ctrl$Drugs !="abine",]
dko_results_ctrl <- read_csv('./Integrated_Drug_Db/DKO_result_CTRL_model.csv')
sko_results_ctrl %>% pivot_longer(!Drugs,names_to = "Functions",values_to = "Growth_Reduction") ->sko_results_ctrl
dko_results_ctrl %>% pivot_longer(!Drugs,names_to = "Functions",values_to = "Growth_Reduction") ->dko_results_ctrl
sko_results_ctrl$Growth_Reduction <- 1-  sko_results_ctrl$Growth_Reduction
dko_results_ctrl$Growth_Reduction <- 1-  dko_results_ctrl$Growth_Reduction
sko_results_ctrl %>% group_by(Drugs) %>% mutate(mean_growth_reduction = mean(Growth_Reduction)) ->sko_results_ctrl
dko_results_ctrl %>% group_by(Drugs) %>% mutate(mean_growth_reduction = mean(Growth_Reduction)) ->dko_results_ctrl
sko_results_ctrl$Drugs <- firstup(sko_results_ctrl$Drugs)
sko_results_ctrl <- sko_results_ctrl[sko_results_ctrl$Functions %in% c('DM_atp_c_','biomass_maintenance'),]
sko_results_ctrl$Functions <- factor(sko_results_ctrl$Functions,c('DM_atp_c_','biomass_maintenance'))

dko_results_ctrl$Drugs <- firstup(dko_results_ctrl$Drugs)
dko_results_ctrl <- dko_results_ctrl[dko_results_ctrl$Functions %in% c('DM_atp_c_','biomass_maintenance'),]
dko_results_ctrl$Functions <- factor(dko_results_ctrl$Functions,c('DM_atp_c_','biomass_maintenance'))
dko_results_ctrl$Drugs <- str_replace(dko_results_ctrl$Drugs, ' ;', " + " )
dko_results_ctrl$Drug1 <- str_split(dko_results_ctrl$Drugs," \\+ ",simplify = T)[,1]
dko_results_ctrl$Drug2 <- firstup(str_split(dko_results_ctrl$Drugs," \\+ ",simplify = T)[,2])
dko_results_ctrl$Drugs <- str_c(dko_results_ctrl$Drug1, " + ",str_to_lower(dko_results_ctrl$Drug2))
dko_results_ctrl$Drugs <- str_replace(dko_results_ctrl$Drugs," \\+ ","\\/")


P_S_ctrl_1 <- ggplot(sko_results_ctrl,aes(y=reorder(Drugs,mean_growth_reduction,decreasing=T),x=Functions))+
  geom_point(aes(size=Growth_Reduction*100),color='steelblue')+
  ylab("Drugs")+  ggtitle("A") +
  theme_classic() +
  scale_size_continuous(range = c(1,4))+
  theme(axis.text.x = element_text(angle = 90,size = 10),
        axis.text.y = element_text(angle = 0,size = 12))+
  guides(size=guide_legend(title = "Growth\nreduction [%]"))
P_S_ctrl_1
P_S_ctrl_2 <- ggplot(dko_results_ctrl,aes(y=reorder(Drugs,mean_growth_reduction,decreasing=T),x=Functions))+
  geom_point(aes(size=Growth_Reduction),color='steelblue')+
  ylab("Combinations")+  ggtitle("B") +
  theme_classic() +
  scale_size_continuous(range = c(1,4))+
  theme(axis.text.x = element_text(angle = 90,size = 10),
        axis.text.y = element_text(angle = 0,size = 12))+guides(size="none")
P_S_ctrl_2
P_S_ctrl_Final <- P_S_ctrl_1 + P_S_ctrl_2+ plot_layout(ncol = 2,tag_level = "new",guides = "collect")#
P_S_ctrl_Final
ggsave(P_S_ctrl_Final,filename="Figure/Figure_S_XX_Control_Model.png", units="in", width=10, height=8,dpi=300,)


### Summarize the xenograft dataset
integrated_xeno <- read_csv("./Drug_DBs_Resources/Integrated_xenografts_data.csv")
bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
drugs_all <- bbb_data[,c('Drugs','Prediction')]
drugs_all$Drugs <- firstup(drugs_all$Drugs)
integrated_xeno$Drugs <- firstup(integrated_xeno$Drugs)
table(unique(integrated_xeno[,c("Drugs","Database")])["Database"])
merged_xeno_drugs <- left_join(drugs_all,integrated_xeno)
table(unique(merged_xeno_drugs[,c("Drugs","Database","Prediction")])[c("Database","Prediction")])
merged_xeno_drugs$Database_dosage <- str_c(merged_xeno_drugs$Database,"\n(",merged_xeno_drugs$Dosage_in_µM, " µM)")

merged_xeno_drugs %>% group_by(Drugs) %>% 
  mutate(Median_Viability_Reduction =median(Viability_Reduction)) %>% 
  group_by(Drugs,Database) %>% 
  mutate(Median_Viability_Reduction_Database =median(Viability_Reduction)) %>% 
  group_by(Drugs,Database_dosage) %>% 
  mutate(Median_Database_dosage =median(Viability_Reduction))->merged_xeno_drugs
n_distinct(merged_xeno_drugs[,c("Drugs")])             
merged_xeno_drugs$Prediction <- factor(merged_xeno_drugs$Prediction,
                                       c("Single",'Combination',"Anti-brain cancer"))
# Define drugs not predicted for GBM

### Add drug targets 
sko_df  <- read_csv('Integrated_Drug_Db/SKO_result_with_Effective_Targets.csv')
sko_df$TYPE2 <- "Single Drug Deletion"
dko_df  = read_csv('Integrated_Drug_Db/DKO_result_with_Effective_Targets.csv')
dko_df$TYPE2 <- "Double Drug Deletion"
dko_df %>% separate_rows(Drugs,sep=' ;')%>%distinct() ->dko_df
sko_df <- rbind(sko_df,dko_df)
sko_df$TYPE2 <- factor(sko_df$TYPE2,c("Single Drug Deletion","Double Drug Deletion"))
sko_df <- sko_df[sko_df$Drugs !="fludarabine",]

gbm_drugs <- firstup(unique(sko_df$Drugs[sko_df$Subtype=='GBM']))
merged_xeno_drugs$Predicted_for_GBM <- "No"
merged_xeno_drugs$Predicted_for_GBM[merged_xeno_drugs$Drugs %in% gbm_drugs] <- "Yes"
merged_xeno_drugs$Predicted_for_GBM[merged_xeno_drugs$Drugs %in% gbm_drugs] <- "Yes"
merged_xeno_drugs$Predicted_for_GBM[merged_xeno_drugs$Prediction %in% "Anti-brain cancer"] <- "Yes"

merged_xeno_drugs$Predicted_for_GBM <- factor(merged_xeno_drugs$Predicted_for_GBM,
                                       c("Yes","No"))
merged_xeno_drugs$Database <- factor(merged_xeno_drugs$Database,c("Bell et al 2018, initial","Bell et al 2018, followup",
                                                                  "Stathias et al 2018"))
merged_xeno_drugs$Database_dosage <- factor(merged_xeno_drugs$Database_dosage,c("Bell et al 2018, initial\n(10 µM)",
                                                                         "Bell et al 2018, followup\n(0.1 µM)",
                                                                         "Bell et al 2018, followup\n(1 µM)",
                                                                         "Bell et al 2018, followup\n(10 µM)",
                                                                  "Stathias et al 2018\n(1 µM)"))
P_5_Xenograft <- ggplot(merged_xeno_drugs[!is.na(merged_xeno_drugs$Viability_Reduction),], 
                   aes(x=Viability_Reduction, #y=reorder_within(Drugs,Median_Viability_Reduction_Database,Database,Dosage_in_µM), 
                       y=reorder_within(Drugs,Median_Database_dosage,Database_dosage), 
                       fill=Prediction,alpha=Predicted_for_GBM,
                       #shape=Subtype,
                       color=Prediction)) + 
  scale_y_reordered()+
  
  geom_violin()+# alpha=0.3
  geom_point()+ #alpha=0.5 #aes(size=Dosage_in_µM)
  geom_point(aes(x=Median_Database_dosage),color='black')+
  #scale_x_continuous(trans='log10',n.breaks = 10)+
  scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  facet_grid(Database_dosage~.,scale="free_y")+
  ylab("Drugs") +
  scale_alpha_discrete(range=c(0.7,0.15))+
  #facet_grid(.~primary_or_metastasis)+
  scale_size_continuous(trans='log10',range=c(1,3))+
  xlab("Growth reduction in GBM xenografts (%)") +theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 10),
        axis.text.y = element_text(angle = 0,size = 11))+
  geom_vline(xintercept=25, linetype="dashed", color = "blue", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  
  ggtitle("A") +
  guides(fill=guide_legend(title='Prediction',ncol = 1,title.position = "top"),
         color=guide_legend(title='Prediction',ncol = 1,title.position = "top"),
         alpha = guide_legend(title='Predicted for GBM',ncol = 1,title.position = "top",
                              override.aes = list(fill = c("forestgreen","lightgreen"))))+
         #size=guide_legend(title='Dosage in µM',ncol = 1,title.position = "top"))+
  #theme(legend.position='bottom')
  theme(legend.background = element_rect(fill="#E1F5FE",size=0.5, linetype="solid",colour ="lightblue"),
        legend.position = c(0.85,0.73))
P_5_Xenograft

### Read xenograft data from literature
xenograft <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S13_inVivo_data")
colnames(xenograft)
xenograft <- xenograft[,c("Drug", "Injected_cell_line","Tested organism","Effect"  ,"Effect with TMZ","Prediction" )]
colnames(xenograft) <-c("Drugs", "Injected_cell_line","Tested_organism"  ,"Effect"  ,"Effect_with_TMZ","Prediction"  )
xenograft$Effect_with_TMZ[is.na(xenograft$Effect_with_TMZ)] <- 'Not tested'
xenograft$Effect[xenograft$Effect %in% c("No effect on tumor growth","No effect on survival")] <-"No effect on tumor growth/survival"
#xenograft$Effect[xenograft$Effect %in% c("Increased survival","Increased radiosensitivity")] <-"Increased radiosensitivity/survival"
xenograft$Effect[xenograft$Effect %in% c("No effect on tumor growth/survival",'Minimal effect on tumor growth')& xenograft$Effect_with_TMZ !='Not tested'] <-xenograft$Effect_with_TMZ[xenograft$Effect %in% c("No effect on tumor growth/survival",'Minimal effect on tumor growth')& xenograft$Effect_with_TMZ !='Not tested']
unique(xenograft$Effect)
xenograft$Effect[xenograft$Effect %in% c("Synergistic","Chemosensitizer")] <- str_c(xenograft$Effect[xenograft$Effect %in% c("Synergistic","Chemosensitizer")],' with TMZ')
unique(xenograft$Effect)

xenograft$Literature <- NA
xenograft$Rank <- 0
xenograft$Rank[str_detect(xenograft$Effect,"Reduce")] <-7
xenograft$Rank[str_detect(xenograft$Effect,"Increased")] <-6

xenograft$Rank[str_detect(xenograft$Effect,"Synergistic")] <-5
xenograft$Rank[str_detect(xenograft$Effect,"Chemosensitizer")] <-4

xenograft$Rank[str_detect(xenograft$Effect,"Minimal")] <-3

xenograft$Rank[str_detect(xenograft$Effect,"Reducing acidosis and hypoxia")] <-2
xenograft$Rank[str_detect(xenograft$Effect,"No effect")] <- 1

xenograft %>% group_by(Drugs)%>% mutate(median_rank = mean(Rank)) -> xenograft
#xenograft$Rank <- 0
xenograft$Prediction <- factor(xenograft$Prediction,c('Single','Combination'))
xenograft$Effect <- factor(xenograft$Effect,c("Reduced tumor growth", "Increased survival"
                                              , "Increased radiosensitivity",
                                              "Synergistic with TMZ" , "Chemosensitizer with TMZ", 
                                              "Minimal effect on tumor growth" ,
                                              "Reducing acidosis and hypoxia" , "No effect on tumor growth/survival"))           


P_5_Xenograft_2 <- ggplot(xenograft,aes(y=reorder(firstup(Drugs),median_rank),x = reorder(Effect,Rank)))+
  #geom_bar(stat = "identity", position=position_dodge(),aes(fill=Source))+
  geom_point(size=3,aes(color=Prediction,fill=Prediction,shape=Tested_organism),
             position = position_dodge(width = 0.8)
  )+
  #geom_text(aes(label=Rank),size=3.5)+
  ggtitle("B") +
  xlab("Drug effect in literature in brain cancer xenografts") +
  ylab("Drugs") +theme_classic() +
  scale_fill_manual(values=c("forestgreen",'#ffa100' ))+
  scale_color_manual(values=c("forestgreen",'#ffa100' ))+
  #facet_wrap(.~Prediction,scales = 'free')+
  #scale_color_manual(values = c("chartreuse4", "green3",'deepskyblue4','deepskyblue1','lightblue','lightred','red'))+
  #scale_fill_manual(values = c("yellow3","violet","white"))+
  #guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(angle = 90,size = 11),
        axis.text.y = element_text(angle = 0,size = 12))+ 
  theme(legend.background = element_rect(fill="#E1F5FE",size=0.5, linetype="solid",colour ="lightblue"),
        legend.position = c(0.25,0.7))

#P5_Final <- P_5_Xenograft + plot_spacer()+ P_5_Xenograft_2+ plot_layout(ncol = 3,widths = c(5, -0.2 ,3))#
P5_Final <- cowplot::plot_grid(P_5_Xenograft, P_5_Xenograft_2, ncol = 2,rel_widths = c(4.75,3.2))
#P5_Final
ggsave(P5_Final,filename="Figure/Figure_5_Xenograft_data.png", units="in", width=12, height=12.2,dpi=300,)

### Read viability reduction in HTS databases and data from literature
Merged_viability <- read_csv("./Merged_Viability_databases.csv")
Merged_viability$Drugs 
colnames(vitro_viability)
colnames(Merged_viability)
intersect(colnames(vitro_viability),colnames(Merged_viability))
vitro_viability <- vitro_viability[,c("Drug","Cell_lines","Reduction_in_Viability" ,"Subtype","Dosage_in_µM")]
vitro_viability$Database <- "Literature"
colnames(vitro_viability) <- colnames(Merged_viability)
vitro_viability$Reduction_in_Viability <- as.numeric(vitro_viability$Reduction_in_Viability)
vitro_viability$Reduction_in_Viability <- vitro_viability$Reduction_in_Viability*100
Merged_viability <- rbind(Merged_viability,vitro_viability)
bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
drugs_all <- bbb_data[,c('Drugs','Prediction')]
drugs_all$Drugs <- firstup(drugs_all$Drugs)
Merged_viability$Drugs <- firstup(Merged_viability$Drugs)

Merged_viability$Drugs[Merged_viability$Drugs =='Fludarabine'] <- "Fludarabine-phosphate"
Merged_viability$Drugs[Merged_viability$Drugs =='Fludarabine phosphate'] <- "Fludarabine-phosphate"
Merged_viability$Drugs[Merged_viability$Drugs =='Doxorubicin hydrochloride'] <- "Doxorubicin"
Merged_viability$Drugs[Merged_viability$Drugs =='5 fluorouracil'] <- "Fluorouracil"
Merged_viability$Drugs[Merged_viability$Drugs =='Doxorubicin hydrochloride'] <- "Doxorubicin"

table(unique(Merged_viability[,c("Drugs","Database")])["Database"])
merged_viability_drugs <- left_join(drugs_all,Merged_viability)
table(unique(merged_viability_drugs[,c("Drugs","Database","Prediction")])[c("Database","Prediction")])
table(unique(Merged_viability[,c('Database','Subtype',"Cell_lines")])[,c('Subtype',"Database")])

merged_viability_drugs %>% group_by(Drugs) %>% 
  mutate(Median_Viability_Reduction =median(Reduction_in_Viability)) %>% 
  group_by(Drugs,Database) %>% 
  mutate(Median_Viability_Reduction_per_database =median(Reduction_in_Viability)) ->merged_viability_drugs

n_distinct(merged_viability_drugs[,c("Drugs")])             
merged_viability_drugs$Prediction <- factor(merged_viability_drugs$Prediction,
                                       c("Single",'Combination',"Anti-brain cancer"))
# Define drugs not predicted for GBM
# gbm_drugs <- firstup(unique(sko_df$Drugs[sko_df$Subtype=='GBM']))
# merged_viability_drugs$Predicted_for_GBM <- "No"
# merged_viability_drugs$Predicted_for_GBM[merged_viability_drugs$Drugs %in% gbm_drugs] <- "Yes"
# merged_viability_drugs$Predicted_for_GBM[merged_viability_drugs$Drugs %in% gbm_drugs] <- "Yes"
# merged_viability_drugs$Predicted_for_GBM[merged_viability_drugs$Prediction %in% "Anti-brain cancer"] <- "Yes"
# 
# merged_viability_drugs$Predicted_for_GBM <- factor(merged_viability_drugs$Predicted_for_GBM,
#                                               c("Yes","No"))
sd(merged_viability_drugs$Dosage_in_µM[!is.na(merged_viability_drugs$Reduction_in_Viability)])
max(merged_viability_drugs$Dosage_in_µM[!is.na(merged_viability_drugs$Reduction_in_Viability)])
merged_viability_drugs$Database <- factor(merged_viability_drugs$Database ,c("PRISM Primary","Literature","Nam et al 2021"))
merged_viability_drugs <- merged_viability_drugs[!is.na(merged_viability_drugs$Reduction_in_Viability),]

P_S_Viability_1 <- ggplot(merged_viability_drugs[merged_viability_drugs$Database=="PRISM Primary",], #[merged_viability_drugs$Database %in% c("PRISM Primary"),]
                        aes(x=Reduction_in_Viability, 
                            #y=reorder_within(Drugs,Median_Viability_Reduction_per_database,Database), 
                            y=reorder(Drugs,Median_Viability_Reduction_per_database), 
                            fill=Prediction,#alpha=Predicted_for_GBM,
                            color=Prediction)) + 
  geom_violin( alpha=0.7)+#
  geom_point(aes(size=Dosage_in_µM),alpha=0.4,position = position_dodge(width = 1))+ #
  geom_point(aes(x=Median_Viability_Reduction_per_database),color='black')+
  #scale_x_continuous(trans='log10',n.breaks = 10)+
  scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  ylab("Drugs") +
  #scale_alpha_discrete(range=c(0.7,0.15))+
  #facet_wrap2(.~Database,scales = "free",ncol = 2,trim_blank = T,shrink = T)+
  #scale_y_reordered()+
  #scale_x_continuous(trans="log10")+
  scale_size_continuous(trans="log10")+#range=c(2,4),breaks = c(2.5,10,500,10000)
  xlab("Viability reduction in the brain cancer cell lines (%)") +theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 11))+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  #ggtitle("Viability reduction in the brain cancer cell lines across two cell viability databases and literature") +
  ggtitle("Primary PRISM") +
  guides(fill=guide_legend(title='Prediction',ncol = 1,title.position = "top"),
         color=guide_legend(title='Prediction',ncol = 1,title.position = "top"),
         alpha = "none",
         size=guide_legend(title='Dosage in µM',ncol = 1,title.position = "top"))+
  #theme(legend.position='bottom') +force_panelsizes(rows = c(5,5,10,0), cols = c(4),respect = T) 
  theme(legend.position = c(0.2,0.2),legend.background = element_blank(),legend.key = element_blank())
P_S_Viability_1

P_S_Viability_2 <- ggplot(merged_viability_drugs[merged_viability_drugs$Database=="Literature",], #[merged_viability_drugs$Database %in% c("PRISM Primary"),]
                          aes(x=Reduction_in_Viability, 
                              #y=reorder_within(Drugs,Median_Viability_Reduction_per_database,Database), 
                              y=reorder(Drugs,Median_Viability_Reduction_per_database), 
                              fill=Prediction,#alpha=Predicted_for_GBM,
                              color=Prediction)) + 
  geom_violin( alpha=0.7)+#
  geom_point(aes(size=Dosage_in_µM),alpha=0.4,position = position_dodge(width = 1))+ #
  geom_point(aes(x=Median_Viability_Reduction_per_database),color='black')+
  #scale_x_continuous(trans='log10',n.breaks = 10)+
  scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  ylab("Drugs") +
  #scale_alpha_discrete(range=c(0.7,0.15))+
  #facet_wrap2(.~Database,scales = "free",ncol = 2,trim_blank = T,shrink = T)+
  #scale_y_reordered()+
  scale_size_continuous(trans="log10")+#range=c(2,4),breaks = c(2.5,10,500,10000)
  xlab("Viability reduction in the brain cancer cell lines (%)") +theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 11))+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  #ggtitle("Viability reduction in the brain cancer cell lines across two cell viability databases and literature") +
  ggtitle("Literature") +
  guides(fill="none",
         color="none",
         alpha = "none",
         size=guide_legend(title='Dosage in µM',ncol = 3,title.position = "top"))+
  #theme(legend.position='bottom') +force_panelsizes(rows = c(5,5,10,0), cols = c(4),respect = T) 
  theme(legend.position = c(0.6,0.15),legend.background = element_blank(),legend.key = element_blank())
P_S_Viability_2

P_S_Viability_3 <- ggplot(merged_viability_drugs[merged_viability_drugs$Database=="Nam et al 2021",], #[merged_viability_drugs$Database %in% c("PRISM Primary"),]
                          aes(x=Reduction_in_Viability, 
                              #y=reorder_within(Drugs,Median_Viability_Reduction_per_database,Database), 
                              y=reorder(Drugs,Median_Viability_Reduction_per_database), 
                              fill=Prediction,#alpha=Predicted_for_GBM,
                              color=Prediction)) + 
  geom_violin( alpha=0.7)+#
  geom_point(aes(size=Dosage_in_µM),alpha=0.4,position = position_dodge(width = 1))+ #
  geom_point(aes(x=Median_Viability_Reduction_per_database),color='black')+
  #scale_x_continuous(trans='log10',n.breaks = 10)+
  scale_fill_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  ylab("Drugs") +
  #scale_alpha_discrete(range=c(0.7,0.15))+
  #facet_wrap2(.~Database,scales = "free",ncol = 2,trim_blank = T,shrink = T)+
  #scale_y_reordered()+
  scale_size_continuous(trans="log10")+#range=c(2,4),breaks = c(2.5,10,500,10000)
  xlab("Viability reduction in the brain cancer cell lines (%)") +theme_classic() +
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 11))+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  #ggtitle("Viability reduction in the brain cancer cell lines across two cell viability databases and literature") +
  ggtitle("Nam et al 2021") +
  guides(fill="none",
         color="none",
         alpha = "none",
         size=guide_legend(title='Dosage in µM',ncol = 1,title.position = "top"))+
  #theme(legend.position='bottom') +force_panelsizes(rows = c(5,5,10,0), cols = c(4),respect = T) 
  theme(legend.position = c(0.6,0.3),legend.background = element_blank(),legend.key = element_blank())

#theme(legend.position = c(0.85,0.62),legend.background = element_blank(),legend.key = element_blank())
design <- matrix(c(1,1,2,3), 2, 2)
P_S_Viability_Final <- cowplot::plot_grid(P_S_Viability_1,(P_S_Viability_2/ P_S_Viability_3),
                                          ncol = 2,rel_heights = c(15,10))
#P_S_Viability_1 <- P_S_Viability_1 + facet_manual(.~Database,scales = "free", design = design)
ggsave(P_S_Viability_Final,filename="Figure/Figure_S_XX_inVitro_viability_reduction.png", units="in", width=11, height=10,dpi=300,)




## Figure 6.b:  Plot the survival data of the clinical trials
trials2 <- read_excel("Supplementary File 2.xlsx",sheet = "Table_S12_Clinical_Trials")
trials2$Arms <- str_replace(trials2$Arms,"Standard therapy","SOC")
colnames(trials2)
colnames(trials2)[colnames(trials2)=="Primary outcome p-value"] <- "P_value"
colnames(trials2)[colnames(trials2)=="Diagnosis abbreviation"] <- "Disease"

# count the numnber of phase1/2 trials in brain cncer
trials2_cnt <- trials2[trials2$Phases %in% c("PHASE2","PHASE3","PHASE1/2") & trials2$`Efficacy as chemotherapy` != "Excluded",]
n_distinct(trials2_cnt$`NCT Number`)

trials2 <- trials2[1:39,]
unique(trials2$Arms)
unique(trials2$Diagnosis)
unique(trials2$Disease)
colnames(trials2)
unique(trials2$`NCT Number`)
unique(trials2$`Intervention Regime`)

trials2$`OS (median in months)`[trials2$`OS (median in months)`=='Not reached'] <- NA
trials2$`OS (median in months)` <- as.double(trials2$`OS (median in months)`)
trials2$`Num evaluable participants`[trials2$`Num evaluable participants`=='No survival data'] <- 2
trials2$`Num evaluable participants` <- as.numeric(trials2$`Num evaluable participants`)

NCTs <- unique(trials2$`NCT Number`)
NCTs <- NCTs[!is.na(NCTs)]
trials2_df <- as.data.frame(trials2)
trials2 <- trials2[rowSums(is.na(trials2)) != ncol(trials2), ]

trials2 %>% replace_na(list(Diagnosis = 'missing', P_value = 'missing',Year=0,
                            Phases = 'missing',`Intervention Regime` ='missing')) -> trials2

for(i in 1:length(NCTs)) {
  nct_i <- NCTs[i]
  diagnosis_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'Diagnosis']))
  disease_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'Disease']))
  
  regime_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'Intervention Regime']))
  phase_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'Phases']))
  pvalue_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'P_value']))
  year_i <- as.vector(na.omit(trials2_df[trials2_df$`NCT Number`==nct_i,'Year']))

  trials2[trials2$`NCT Number`==nct_i,'Diagnosis'] <- diagnosis_i
  trials2[trials2$`NCT Number`==nct_i,'Disease'] <- disease_i
  
  trials2[trials2$`NCT Number`==nct_i,'Intervention Regime'] <- regime_i
  trials2[trials2$`NCT Number`==nct_i,'Phases'] <- phase_i
  #trials2[trials2$`NCT Number`==nct_i,'P_value'] <- as.list(pvalue_i)
  #trials2[trials2$`NCT Number`==nct_i,'Year'] <- year_i
}


trials2 %>%  pivot_longer(cols =c("OS (median in months)","PFS (median in months)") ,
                          names_to = "Measure",values_to = "Measure_value") ->trials2

trials2 %>%group_by(`NCT Number`) %>% mutate(mean_survival = mean(Measure_value,na.rm=T)) ->trials2

trials2$NCT_Year <- str_c(str_split(trials2$`NCT Number`,"_",simplify = T)[,1],"\n(",trials2$Year,")")
unique(trials2$`NCT Number`)
unique(trials2$NCT_Year)

trials2$Diagnosis
failed_trials <- c("Cladribine_1999","Melphalan_1988","Gemcitabine_2000")
trials2[trials2$`NCT Number` %in% failed_trials,'Num evaluable participants'] <- trials2[trials2$`NCT Number` %in% failed_trials,'Number of participants']
trials2[trials2$`NCT Number` %in% failed_trials,'Arms'] <- "Non-effective"
trials2[trials2$`NCT Number` %in% failed_trials,'Measure_value'] <- 0
trials2$`Intervention Regime` <- str_replace(trials2$`Intervention Regime`,"After 1st line ","")
trials2$`Intervention Regime` <- str_replace(trials2$`Intervention Regime`,"\\)","")
trials2$`Intervention Regime` <- str_replace(trials2$`Intervention Regime`,"\\(","")

unique(trials2$`Intervention Regime`)
unique(trials2$Diagnosis)
trials2 %>% group_by(`NCT Number`) %>%
  mutate(n_Arms = n_distinct(Arms)) -> trials2
#trials2$`Intervention Regime` <- factor(trials2$`Intervention Regime` ,c("Monotherapy","TMZ / RT","RT","TMZ" ,  "TMZ + RT","bevacizumab", "Immunotherapy"))
#trials2 <- as_tibble(trials2)
#trials2 <- trials2[! is.na(trials2$Measure_value),]

trials2 <- trials2[trials2$`Intervention Regime`!='Immunotherapy',]
#trials2 <- trials2[!is.na(trials2$Measure_value),]

trials2$Diagnosis[trials2$Diagnosis=="Neuroblastoma and other childhood solid tumour"] <-"Neuroblastoma"
trials2$Diagnosis[trials2$Diagnosis=="Melanoma Brain Metastasis"] <-"Melanoma BM"
trials2$`Intervention Regime`[trials2$`Intervention Regime`=="Vincristine + Lomustine + Procarbazine"] <-"PCV"

unique(trials2$`Intervention Regime`)
unique(trials2$Disease)
trials2$Measure[trials2$Measure == "PFS (median in months)" ] <- "Median PFS" 
trials2$Measure[trials2$Measure == "OS (median in months)" ] <- "Median OS" 

trials2$Disease <- factor(trials2$Disease ,
                          c("rGBM","nGBM","rLGG" ,"rGlioma","nGlioma","DIPG","Neuroblastoma",
                            "Meningioma"))#,"Melanoma BM"

PS_XX_trials <- ggplot(trials2,
                       aes(y=reorder(NCT_Year,mean_survival),
                           label=Arms,
                 shape=Phases,
                 color=`Intervention Regime`,
                 #color=`Intervention Regime`,
                 x =Measure_value))+
  geom_point(aes(size=`Num evaluable participants`),position = position_dodge2(width = 1,reverse = T))+
  geom_text_repel(position = position_dodge2(width = 2,reverse = T),size=5)+
  #geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  #geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  facet_grid(Disease~Measure,scales = 'free_y')+
  ggtitle("") +
  #xlab("") +
  force_panelsizes(rows = c(7,3,2,3,5,1,2,2),
                   cols = c(3, 3)) +
  xlab("Survival") +
  ylab("Clinical trials") +
  theme_classic() +
  scale_x_continuous(trans = 'log10')+
  scale_size_continuous(trans = 'log10',range = c(3,7))+
  #scale_x_continuous(trans = 'log10')+
  #scale_size_continuous(trans = 'log10',range = c(3,7))+
  #theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))+
   scale_shape_manual(values = seq(1,13,1)) 

PS_XX_trials


ggsave(PS_XX_trials,filename="Figure/FigS_XX_Trials.png", units="in", width=12.5, height=14, dpi=300)


##Only two-arms trials
trials_effective_drugs <- c("Celecoxib","Eflornithine","Valganciclovir","Fotemustine")
trials2_effective_ncts <-trials2$`NCT Number`[str_detect(trials2$Arms,paste0(trials_effective_drugs,collapse = "|"))]

trials2_effective <-trials2[trials2$`NCT Number` %in% trials2_effective_ncts,]
trials2_effective <- trials2_effective[trials2_effective$Disease %in% c("rGBM","nGBM","rLGG","LGG","AG"),]
trials2_effective <- trials2_effective[trials2_effective$n_Arms>1,]
trials2_effective$Prediction <- "NA"


bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
drugs_all <- bbb_data[,c('Drugs','Prediction')]

drugs_single <- firstup(drugs_all$Drugs[drugs_all$Prediction=="Single"])
drugs_comb<- firstup(drugs_all$Drugs[drugs_all$Prediction=="Combination"])
drugs_approved <- c("SOC","TMZ","Bevacizumab","PCV")
trials2_effective$Prediction[str_detect(trials2_effective$Arms,paste0(drugs_approved,collapse = "|"))] <- "Anti-brain cancer"
trials2_effective$Prediction[str_detect(trials2_effective$Arms,paste0(drugs_single,collapse = "|"))] <- "Single"
trials2_effective$Prediction[str_detect(trials2_effective$Arms,paste0(drugs_comb,collapse = "|"))] <- "Combination"
trials2_effective$Prediction <- factor(trials2_effective$Prediction,
                                            c("Single",'Combination',"Anti-brain cancer"))
trials2_effective$Arms_pvalue <- str_c(trials2_effective$Arms,"\n(p = ",trials2_effective$P_value,")")
trials2_effective$Arms_pvalue[trials2_effective$P_value=="missing"] <- trials2_effective$Arms[trials2_effective$P_value=="missing"]



P6_2 <- ggplot(trials2_effective,
               aes(y=reorder(NCT_Year,mean_survival),
                   label=Arms_pvalue,
                   shape=Phases,
                   color= Prediction,
                   #color=`Intervention Regime`,
                   #color=`Intervention Regime`,
                   x =Measure_value))+
  geom_point(aes(size=`Num evaluable participants`),position = position_dodge2(width = 1,reverse = T))+
  geom_text_repel(position = position_dodge2(width = 3,reverse = T,padding = 0.05),size=4.5,force_pull = 1)+
  #geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  #geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  facet_grid(Disease~Measure,scales = 'free')+
  ggtitle("") +
  xlab("Survival (in months)") +
  ylab("Clinical trials") +
  scale_color_manual(values=c("forestgreen",'#ffa100','steelblue'))+
  force_panelsizes(rows = c(5,3.5,2),
                   cols = c(4.5, 3)) +
  theme_classic() +
  ggtitle("B")+
  
  scale_x_continuous(trans = 'log10')+
  scale_size_continuous(trans = 'log10',range = c(3,7),guide = guide_legend(title = "Number of\nevaluable patients"))+
  #theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 0,size = 14),
        axis.text.y = element_text(angle = 0,size = 14),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12)
  )+
  scale_shape_manual(values = seq(1,13,1)) 

P6_2

#P6_final <- cowplot::plot_grid(P6_1, P6_2, ncol = 1,rel_heights = c(5,9))

#ggsave(P6_final,filename="Figure/Figure_6_trials.png", units="in", width=11, height=11, dpi=300)


# Figure 6.a: The number of effective, ineffective, untested drugs in vitro, in vivo, and in clinical trials
effective_single <-  read_excel("Supplementary File 2.xlsx",sheet = "Table_S9_SingleD_Selection")
colnames(effective_single)
effective_comb <- read_excel("Supplementary File 2.xlsx",sheet = 3)
colnames(effective_comb)
selected_cols <- c("Drugs" ,"Target pathway","Clinical trials as combination in 2-arms, phase 1/2 or higher",
                   "Overall in/ex vivo xenografts results","Overall in vitro results"  )
effective_single <- effective_single[1:33,c(1,2,5,6,7)]
effective_comb <- effective_comb[1:19,selected_cols]
effective_single$Prediction <- "Single"
effective_comb$Prediction <- "Combination"
newcolnames <- c("Drugs" ,"Indication","Clinical_trials","Xenografts","Invitro" ,"Prediction" )
colnames(effective_single) <- newcolnames
colnames(effective_comb) <- newcolnames
effective_df <- rbind(effective_single,effective_comb)
effective_df_keep <- effective_df
effective_df$Indication[effective_df$Indication!= "Anticancer"] <- "Non anticancer"
effective_df$Indication[effective_df$Drugs== "resveratrol"] <- "Non anticancer"
effective_df$Drugs <- firstup(effective_df$Drugs)
effective_df$Clinical_trials[effective_df$Drugs %in% c("Celecoxib","Eflornithine","Valganciclovir","Fotemustine")] <- "Effective"
effective_df$Clinical_trials[effective_df$Drugs %in% c("Gemcitabine","Melphalan","Cladribine","Mercaptopurine","Fluorouracil")] <- "Ineffective"
effective_df$Clinical_trials[!effective_df$Clinical_trials %in% c("Effective","Ineffective","Untested")] <- "Untested"
effective_df %>% pivot_longer(cols =c("Clinical_trials","Xenografts","Invitro") ,
                              names_to = "Evidence",values_to = "Results") -> effective_df_longer
effective_df_longer %>% group_by(Prediction,Indication,Evidence,Results) %>%
  summarise(N_drugs_per_indication = n_distinct(Drugs))%>% 
  group_by(Prediction,Evidence,Results) %>%
           mutate( N_drugs_total = sum(N_drugs_per_indication)) -> effective_df_summ
# Filling zero summary
colnames(effective_df_summ)
effective_df_summ_2 <- data.frame(Prediction=c("Single","Combination"),
                                  Indication = c("Non anticancer","Non anticancer"),
                                  Evidence  = c("Invitro","Clinical_trials")  ,
                                  Results  = c("Untested","Ineffective") ,
                                  N_drugs_per_indication= c(0,0) ,
                                  N_drugs_total = c(0,0)
                                    )
effective_df_summ <- rbind(effective_df_summ,effective_df_summ_2)
effective_df_summ$Evidence_ <- effective_df_summ$Evidence
effective_df_summ$Evidence_[effective_df_summ$Evidence_ =="Invitro"]<- "In vitro"
effective_df_summ$Evidence_[effective_df_summ$Evidence_ =="Xenografts"]<- "Xenografts"
effective_df_summ$Evidence_[effective_df_summ$Evidence_ =="Clinical_trials"]<- "Clinical trials (Phase 2)"


effective_df_summ$Prediction <- factor(effective_df_summ$Prediction,c("Single","Combination"))
effective_df_summ$Evidence <- factor(effective_df_summ$Evidence,c("Invitro","Xenografts","Clinical_trials"))
effective_df_summ$Evidence_ <- factor(effective_df_summ$Evidence_,
                                      c("In vitro","Xenografts","Clinical trials (Phase 2)"))

effective_df_summ$Results <- factor(effective_df_summ$Results,rev(c("Effective","Ineffective","Untested")))
effective_df_summ$Indication <- factor(effective_df_summ$Indication,c("Non anticancer","Anticancer"))
effective_df_summ <- effective_df_summ[! (effective_df_summ$Prediction == "Combination" & effective_df_summ$N_drugs_total==0), ]
P6_1<- ggplot(effective_df_summ,aes(y=Results,
                               fill=Results,
                               alpha=Indication,
                               x = N_drugs_per_indication))+
  geom_bar(stat = "identity")+
  #geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  #geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  xlab("Number of drugs") +
  ggtitle("A")+
  ylab("Evidence results") +theme_classic() +
  #scale_x_continuous(limits=c(0,110))+
  facet_grid(Prediction~Evidence_)+
  #scale_fill_brewer(palette="Dark2")+
  geom_text(aes(x=N_drugs_total,label=N_drugs_total,color=Results),size=4,nudge_x = 3,alpha=1)+
  
  scale_fill_manual(values = c("#505050",'red', "forestgreen"))+#"#105d46"
  scale_color_manual(values = c("#505050",'red', "forestgreen"))+#"#105d46"
  
  scale_alpha_discrete(range=c(0.3,1),guide = guide_legend(title ="Approved indication"))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill='none',color='none')+
  theme(axis.text.x = element_text(angle = 0,size = 14),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 14),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12)
  )
  #theme(legend.position = c(0.55,0.6),legend.background = element_blank(),legend.key = element_blank())
P6_1
P6_final <- cowplot::plot_grid(P6_1, P6_2, ncol = 1,rel_heights = c(5,9))

ggsave(P6_final,filename="Figure/Figure_6_trials.png", units="in", width=11, height=11, dpi=300)
ggsave(P6_final,filename="Figure/Figure_6_trials.pdf", units="in", width=11, height=11, dpi=300)


#P_S_Viability_1+  P_S_Viability_2 + plot_layout(ncol = 2,guides = "collect") +scale_size_continuous(trans="log10")#+
#ggarrange(P_S_Viability_1, P_S_Viability_2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
#plot_grid(P_S_Viability_1,  P_S_Viability_2,)
#ggsave(P_S_Viability_1,filename="Figure/Figure_S_XX_inVitro_viability.pdf", units="in", width=12, height=11,dpi=300,)


# sko_df_ <- distinct(sko_df_longer[,c('Drugs','TYPE2')])
# sko_df_ <- left_join(sko_df_,xenograft)
# bell_df <- left_join(sko_df_,Ranking_df_3)
# 
# bell_df <- left_join(bell_df,Ranking_df_4[,c('Drugs','Bell_etal_Folllowup')])
# #bell_df <- left_join(bell_df,xenograft)
# 
# bell_df %>% pivot_longer(cols =c("Bell_etal_Initial","Bell_etal_Folllowup","Literature") ,
#                          names_to = "Source",values_to = "Rank")  -> bell_df
# 
# 
# bell_df %>% group_by(Drugs) %>% mutate(median_rank = median(Rank,na.rm = T)) -> bell_df
# bell_df[bell_df$Source=='Literature','Rank'] <- 0
# bell_df[bell_df$Source=='Literature','median_rank'] <- 0
# bell_df[bell_df$Source!='Literature','Effect'] <- NA
# bell_df[bell_df$Source!='Literature','Xenograft_or_Animal'] <- NA
# 
# 
# min(bell_df$median_rank)
# bell_df$median_rank[bell_df$Drugs %in% c("fotemustine","azathioprine",'eflornithine')] <- -1
# bell_df$median_rank[bell_df$Drugs %in% c("cannabidiol")] <- -0.5
# bell_df$median_rank[bell_df$Drugs %in% c("topiramate")] <- 150
# bell_df$median_rank[bell_df$Drugs %in% c("inositol","hydroxyurea",'acetazolamide')] <- 200
# #bell_df$Effect <- as.character(bell_df$Effect)
# #bell_df$Effect[is.na(bell_df$median_rank)] <- "Not tested"
# #bell_df$Rank[is.na(bell_df$median_rank)] <- 0
# 
# bell_df$median_rank[is.na(bell_df$median_rank)] <- 300
# 
# bell_df$Effect <- factor(bell_df$Effect,c("Reduced tumor growth", "Minimal effect on tumor growth", "Reducing acidosis and hypoxia","No effect on tumor growth/survival",NA   ))
# 
# #bell_df$median_rank[bell_df$Effect %in% c("No effect on tumor growth","No effect on survival")] <-"No effect on tumor growth/survival"
# 
# colnames(bell_df)
# colnames(xenograft)
# ggplot(bell_df,aes(y=fct_reorder(Drugs,desc(median_rank)),x = Rank))+
#   geom_bar(stat = "identity", position=position_dodge(),aes(fill=Source))+
#   geom_point(aes(shape=Xenograft_or_Animal,color=Effect,fill=Source),size=3,position = position_dodge(width = 0.8))+
#   #geom_text(aes(label=Rank),size=3.5)+
#   ggtitle("Xenograft data for the predicted drugs from literature and drug screening") +
#   xlab("Rank in nth percentile among all screened drugs") +
#   ylab("Drugs") +theme_classic() +
#   facet_wrap(.~TYPE2,scales = 'free')+
#   scale_color_manual(values = c("chartreuse4", "green3",'deepskyblue4','red','grey'))+
#   scale_fill_manual(values = c("yellow3","violet","white"))+
#   guides(color=guide_legend(title="Effect alone in Literature"))+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12))



## Rank essential genes and drug targets by DepMap dependency

depmap_dep = read_csv('./CRISPR/DepMap_22Q1/CRISPR_gene_dependency.csv')
depmap_cellinfo <- read.csv("CRISPR/DepMap_22Q1/sample_info_braincancer.csv")

depmap_dep <- as.data.frame(depmap_dep)
rownames(depmap_dep) <- depmap_dep$DepMap_ID

depmap_dep <- depmap_dep[depmap_dep$DepMap_ID %in% depmap_cellinfo$DepMap_ID,]

depmap_dep_longer <- depmap_dep %>% pivot_longer(!DepMap_ID,names_to = "Targets",values_to = "Dependancy_Probability")
depmap_dep_longer$Targets <-  str_split(depmap_dep_longer$Targets," ",simplify = T)[,1]

#metabolic_genes <-  read_csv('./Generic_Models/Recon3D_HMR2_Human1_genes.csv')
#metabolic_genes <- metabolic_genes$SYMBOL
#depmap_scr_dep_gbm <- depmap_scr_dep_gbm[depmap_scr_dep_gbm$SYMBOL %in%metabolic_genes, ]
#length(unique(depmap_scr_dep_gbm$DepMap_ID))
depmap_dep_longer <- left_join(depmap_dep_longer,depmap_cellinfo)
unique(depmap_dep_longer$primary_or_metastasis)

depmap_dep_longer %>% 
  group_by(Subtype) %>%
  mutate(n_cells= n_distinct(DepMap_ID), # number of cell lines by subtype
         Subtype_n_cells =  str_c(Subtype," (n=", n_cells,")"))  %>%
  group_by(Targets,Subtype,Subtype_n_cells) %>%# 
  summarise(median_prob=median(Dependancy_Probability,na.rm = T))%>%
  dplyr::arrange(median_prob)->depmap_dep_longer_rank

depmap_dep_longer_rank %>% group_by(Subtype_n_cells) %>% 
  mutate(rank =  ntile(-median_prob,n_distinct(Targets))) -> depmap_dep_longer_rank
depmap_dep_longer_rank <- depmap_dep_longer_rank[!is.na(depmap_dep_longer_rank$Subtype),]
depmap_dep_longer_rank_placeholder <- depmap_dep_longer_rank

## Merge with the drug targets
colnames(sko_placeholder__)  
sko_placeholder <- sko_placeholder__[,c("SYMBOL","Drugs", "Subtype" , "TYPE2" , "Predicted_Essential")]
unique(sko_placeholder$Predicted_Essential)
unique(sko_placeholder$TYPE2)
#sko_placeholder$TYPE2 <- levels(sko_placeholder$TYPE2)

sko_placeholder$Predicted_Essential[sko_placeholder$Predicted_Essential =="Non-essential drug targets"] <- "No"
sko_placeholder$Predicted_Essential[sko_placeholder$Predicted_Essential =="Essential drug targets"] <- "Yes"
#sko_placeholder$TYPE2[sko_placeholder$TYPE2 =="Single Drug Deletion"] <- "Single drug target"
#sko_placeholder$TYPE2[sko_placeholder$TYPE2 =="Double Drug Deletion"] <- "Combination drug target"
levels(sko_placeholder$TYPE2) <- c("Single drug target","Combination drug target")

essential_df <- data_longer_[data_longer_$Predicted_Essential==1,c("Symbol","Subtype")] #data_longer_[data_longer_$Predicted_Essential==1,]
colnames(essential_df)[1] <- c("SYMBOL")
essential_df$Drugs <- "_"
essential_df$TYPE2 <- "Essential gene"
essential_df$Predicted_Essential <- "Yes"
colnames(sko_placeholder)  

sko_placeholder_ <- rbind(sko_placeholder,essential_df[colnames(sko_placeholder)]) # merge the essential genes and drug targets
n_distinct(sko_placeholder_$SYMBOL)


integrated_db_all <- read_csv("./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv")
integrated_db_all$Drugs <- firstup(integrated_db_all$name)

bbb_data <- read_excel("Supplementary File 2.xlsx", sheet = "Table_S16_CSF_Bioavailability")
approved_df <- bbb_data[bbb_data$Prediction=="Anti-brain cancer",c("Drugs")]
approved_df$Drugs <- firstup(approved_df$Drugs)
approved_df <- unique(left_join(approved_df,integrated_db_all[,c("Drugs","SYMBOL")]))
approved_df$TYPE2 <- "Anti-brain cancer target"
approved_df$Predicted_Essential <- "No"
approved_df$Subtype <- "AST:GBM:ODG"
approved_df %>% separate_rows(Subtype,sep=':')%>%distinct()->approved_df

n_distinct(approved_df$SYMBOL)
intersect(approved_df$SYMBOL,essential_df$SYMBOL)
intersect(approved_df$SYMBOL,sko_placeholder$SYMBOL)

sko_all <- rbind(sko_placeholder_,approved_df[colnames(sko_placeholder)])
sko_all <- unique(sko_all)
colnames(sko_all)[1] <- "Targets"
unique(sko_all$TYPE2)

colnames(depmap_dep_longer_rank)
all_genes <- unique(sko_all$Targets)

# merging on the ranking to know the rank on non predicted genes 
drugTypes <- c("Single drug target","Combination drug target","Essential gene",  "Anti-brain cancer target")

depmap_dep_longer_rank$TYPE2 <- paste0(drugTypes,collapse = ":")
depmap_dep_longer_rank %>% separate_rows(TYPE2,sep=':')%>%distinct()->depmap_dep_longer_rank

depmap_dep_longer_rank <- left_join(depmap_dep_longer_rank,sko_all) 
# A placeholder for ranked depmap genes to evaluate other gene sets

depmap_dep_longer_rank$TYPE2 <- factor(depmap_dep_longer_rank$TYPE2,drugTypes)
                                          
#Keeping only genes in the previous four classes
depmap_dep_longer_rank <- depmap_dep_longer_rank[depmap_dep_longer_rank$Targets %in% all_genes,]
depmap_dep_longer_rank <- depmap_dep_longer_rank[depmap_dep_longer_rank$Subtype != "Non-glioma",]
depmap_dep_longer_rank <- depmap_dep_longer_rank[!(is.na(depmap_dep_longer_rank$Predicted_Essential) & depmap_dep_longer_rank$TYPE2 =="Anti-brain cancer target"), ]
depmap_dep_longer_rank %>% group_by(Targets) %>%
  mutate(to_drop = TYPE2 %in% TYPE2[!is.na(Predicted_Essential)]) -> depmap_dep_longer_rank
depmap_dep_longer_rank <- depmap_dep_longer_rank[depmap_dep_longer_rank$to_drop==T,]
  
## Approved_or_Predicted
# Approved 
# Predicted

## Prediction
# Anti-brain cancer target
# Drug target
# Essential gene


depmap_dep_longer_rank %>% group_by(Targets,Subtype) %>%# 
  mutate(rank=max(rank),median_prob=max(median_prob))%>% distinct() ->depmap_dep_longer_rank
depmap_dep_longer_rank$rank_1 <- 0
depmap_dep_longer_rank <- depmap_dep_longer_rank[!is.na(depmap_dep_longer_rank$median_prob),]
depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$Subtype=='GBM'] <- depmap_dep_longer_rank$rank[depmap_dep_longer_rank$Subtype=='GBM']
#depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$Subtype!='GBM'] <- 20000
#depmap_dep_longer_rank$Tested_InVitro[depmap_dep_longer_rank$Tested_InVitro %in% c("Selected","Not selected")] <- "Predicted drug target"
n_distinct(depmap_dep_longer_rank$Targets)
n_distinct(depmap_dep_longer_rank$Drugs)


#depmap_dep_longer_rank$Gene_type <- factor(depmap_dep_longer_rank$Gene_type,
#                                           c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
#                                             "Anti-melanoma target","NO-related genes"))

depmap_drugs<- depmap_dep_longer_rank
#depmap_dep_longer_rank <- distinct(depmap_dep_longer_rank[,!colnames(depmap_dep_longer_rank) %in% c("Drugs","Tested_InVitro")])

depmap_dep_longer_rank <- depmap_dep_longer_rank%>% group_by(Targets,TYPE2) %>% 
  mutate(median_prob_gbm = unique(median_prob[Subtype=='GBM']))

depmap_dep_longer_rank$Predicted_per_subtype <- "No"
depmap_dep_longer_rank$Predicted_per_subtype[!is.na(depmap_dep_longer_rank$Predicted_Essential)] <- "Yes"


shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

# Make font map for the essential genes
depmap_essential_bold <- depmap_dep_longer_rank[depmap_dep_longer_rank$Subtype=='GBM',]  |> 
  select(Targets,rank_1,Predicted_Essential)|>
  distinct() |>
  dplyr::arrange(desc(rank_1))
depmap_essential_bold$Predicted_Essential[depmap_essential_bold$Predicted_Essential=='Yes'] <- "bold"
depmap_essential_bold$Predicted_Essential[depmap_essential_bold$Predicted_Essential=='No'] <- "plain"
depmap_dep_longer_rank$Subtype_n_cells <- factor(depmap_dep_longer_rank$Subtype_n_cells ,c("GBM (n=44)","AST (n=20)","ODG (n=2)"))
variant_genes <- unique(depmap_dep_longer_rank$Targets[depmap_dep_longer_rank$Predicted_per_subtype=="No"])
# 
# # Removing genes not targetedd with OR rule
# sko_df_longer_not_targeted %>% separate_rows(Rule,sep = " \\| ")-> sko_df_longer_not_targeted_
# sko_df_longer_not_targeted_$Rule <- str_replace(sko_df_longer_not_targeted_$Rule,"\\(","")
# sko_df_longer_not_targeted_$Rule <- str_replace(sko_df_longer_not_targeted_$Rule,"\\)","")
# sko_df_longer_not_targeted_$Drugs <- firstup(sko_df_longer_not_targeted_$Drugs)
# unique(sko_df_longer_not_targeted_$Rule)
# depmap_dep_longer_rank$Targets
# depmap_dep_longer_rank <- left_join(depmap_dep_longer_rank,sko_df_longer_not_targeted_[,c("Drugs","Subtype","Rule","Rule_targeted")],
#                                     by=c("Targets"="Rule","Drugs"="Drugs","Subtype"="Subtype"))
# depmap_dep_longer_rank <- depmap_dep_longer_rank[is.na(depmap_dep_longer_rank$Rule_targeted),]

depmap_dep_longer_rank[,colnames(depmap_dep_longer_rank)!= "TYPE2"] %>% 
  summarise(Targets,median_prob_gbm) %>%
  distinct( )%>%
  dplyr::arrange(-median_prob_gbm)->depmap_dep_longer_genes
  
gene_len <- n_distinct(depmap_dep_longer_genes$Targets) 
genes_sorted_1 = depmap_dep_longer_genes$Targets[1:(gene_len%/%2) ]
genes_sorted_2 = depmap_dep_longer_genes$Targets[((gene_len%/%2)+1) :gene_len]

P_XX_depmap_1 <- ggplot(depmap_dep_longer_rank[depmap_dep_longer_rank$Targets %in%  genes_sorted_1,],
                      aes(y=reorder(Targets,median_prob_gbm),
                                           fill=TYPE2,
                                           alpha=Predicted_per_subtype,
                                           x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  #ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the glioma cell lines out of 17386 genes") +
  ggtitle("A") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug targets/ essential genes") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue' ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue'))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 8)
  )+
  theme(legend.position = c(0.9,0.35),legend.background = element_blank(),legend.key = element_blank())
P_XX_depmap_1

P_XX_depmap_2 <- ggplot(depmap_dep_longer_rank[depmap_dep_longer_rank$Targets %in%  genes_sorted_2,],
                        aes(y=reorder(Targets,median_prob_gbm),
                            fill=TYPE2,
                            alpha=Predicted_per_subtype,
                            x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 20)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("B") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug targets/ essential genes") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue' ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue'))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill='none',color='none',alpha="none")+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 8.5)
  )
P_XX_depmap_2
P_XX_depmap_final <- cowplot::plot_grid(P_XX_depmap_1, P_XX_depmap_2, ncol = 2,rel_widths = c(5,2.5))
P_XX_depmap_final
ggsave(P_XX_depmap_final,filename="Figure/FigS_XX_DepMap_Rank.png", units="in", width=13, height=11.8, dpi=300)

P_XX_depmap_diff <- ggplot(depmap_dep_longer_rank[depmap_dep_longer_rank$Targets %in% variant_genes,],
                      aes(y=reorder(Targets,desc(rank_1)),
                                                    fill=TYPE2,
                                                    alpha=Predicted_per_subtype,
                                                    x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the glioma cell lines out of 17386 genes") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug Targets") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue' ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue'))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 9)
  )+
  theme(legend.position = c(0.55,0.6),legend.background = element_blank(),legend.key = element_blank())
P_XX_depmap_diff
ggsave(P_XX_depmap_diff,filename="Figure/FigS_XX_DepMap_Rank_variable.png", units="in", width=8, height=6, dpi=300)

## Genes shared beween ABCs and predicted
shared_genes <- intersect(approved_df$SYMBOL,sko_placeholder$SYMBOL)

depmap_dep_longer_rank_shared <- depmap_dep_longer_rank[depmap_dep_longer_rank$Targets %in% shared_genes,]

P_XX_depmap_shared <- ggplot(depmap_dep_longer_rank_shared,
                           aes(y=reorder(Targets,desc(rank_1)),
                               fill=TYPE2,
                               alpha=Predicted_per_subtype,
                               x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the glioma cell lines out of 17386 genes") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug Targets") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue' ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue'))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 9)
  )+
  theme(legend.position = c(0.55,0.6),legend.background = element_blank(),legend.key = element_blank())
P_XX_depmap_shared
ggsave(P_XX_depmap_shared,filename="Figure/FigS_XX_DepMap_Rank_shared.png", units="in", width=8, height=6, dpi=300)



## Check the FVA for all reactions
fva_all <- read_csv('./Sample_models/Subtypes_Models_FVA_All_Reactions.csv') # Fluxes of all reactions

fva_all %>% pivot_longer(cols = c( "min_AST", "max_AST", "min_GBM", "max_GBM" ,"min_ODG" ,"max_ODG"),
                         names_to = 'Model',values_to = 'Flux') -> fva_all
fva_all %>% group_by(Rxn) %>%
  mutate(to_keep = if_else(is.na(Flux),'remove','keep')) -> fva_all # | abs(Flux)<1

fva_all %>% group_by(Rxn) %>% filter(any(to_keep=='keep')) -> fva_all

#fva_all <- na.omit(fva_all)
#fva_all <- fva_all[fva_all$Flux>1e-4,]
unique(fva_all$Rxn)
fva_all$Subtype <- str_to_upper(str_split(fva_all$Model,'_',simplify = T)[,2])
fva_all$FluxType <- str_split(fva_all$Model,'_',simplify = T)[,1]
fva_all %>% group_by(Rxn) %>%
  mutate(average_flux = mean(abs(Flux), na.rm=TRUE)) -> fva_all
# Ceiling the measure to the 2nd decimal
fva_all$Flux_ <- lapply(as.numeric(fva_all$Flux) , function(x) signif.ceiling(x, 3))
fva_all$newname <- str_c(fva_all$Subtype,'_',fva_all$FluxType)
fva_all$FluxType[fva_all$FluxType=='min'] <- 'Min'
fva_all$FluxType[fva_all$FluxType=='max'] <- 'Max'
#fva_all$FluxType <- factor(fva_all$FluxType ,levels = c('Min','Max'))

#va_all$Formula <- str_replace(fva_all$Formula,'  <=>','')
#fva_all %>% group_by(Formula, Subtype) %>%
#  mutate(Reaction=ifelse(sum(Flux)==0,'Inactive','Release')) ->fva_all
fva_all$Reaction <- NA
#fva_all$Reaction[fva_all$ymin=='NaN' & fva_all$ymax=='NaN' ]  <- 'Inactive'
#fva_all$Reaction[fva_all$ymin==0 & fva_all$ymax==0]  <- 'Inactive'

#fva_all$Reaction[fva_all$ymin=>-1e-12 & fva_all$ymax<=1e-12 ]  <- 'Inactive'

#fva_all$Reaction[fva_all$ymin=>0 & fva_all$ymax>0 ]  <-'Release' 3
#fva_all$Reaction[fva_all$ymax<=0 & fva_all$ymin<0 ]  <- 'Uptake'
#fva_all$Reaction[fva_all$ymin<0 & fva_all$ymax>0 ]  <-'Uptake/Release' 
#fva_all$Reaction[fva_all$ymin=>-1e-12 & fva_all$ymax<=1e-12 ]  <- 'Inactive'

fva_all$Reaction[fva_all$Flux==0 | fva_all$Flux =='NaN']  <- 'Inactive'
#fva_all$Reaction[fva_all$ymin=>0 & fva_all$ymax>0 ]  <-'Release' 3
fva_all$Reaction[fva_all$Flux <0 ]  <- 'Uptake'
fva_all$Reaction[fva_all$Flux >0 ]  <-'Release' 
fva_all$VMH_ID <- str_split(fva_all$Rxn,'EX_',simplify = T)[,2]
fva_all$VMH_ID <- str_split(fva_all$VMH_ID,'\\[',simplify = T)[,1]
fva_all$VMH_ID <- str_to_upper(fva_all$VMH_ID)


fva_all <- left_join(fva_all,class_df)
fva_all$super_class[is.na(fva_all$super_class)] <- 'Undefined'
unique(fva_all$Formula[is.na(fva_all$super_class)])

fva_all %>% group_by(Formula) %>% 
  mutate(rank = sum(abs(Flux),na.rm = T))%>% 
  
  group_by(class) %>% 
  mutate(rank_2 = sum(abs(Flux),na.rm = T))-> fva_all #

# rank for narrow-bound reactions with differences between the three sub types
# 1. calculate the difference between min and max in each subtype
# 2.  multiply the difference by the sum(abs(flux)) in each subtype
# 3. take the difference of 2 of GBM - (AST + ODG)
fva_all %>%  group_by(Rxn, Subtype)   %>% 
  mutate(rank_3_1 = Flux[FluxType=="Max"] -  Flux[FluxType=="Min"])  %>% # the smaller, the higher the difference
  group_by(Rxn, Subtype) %>% 
  mutate(rank_3_2 =  sum(abs(Flux))/unique(rank_3_1) )  %>% 
  group_by(Rxn) %>% 
  mutate(rank_3 =  sum(rank_3_2[Subtype !="GBM"] )/unique(rank_3_2[Subtype =="GBM"]))-> fva_all 

## Adding the GPR rules
recon_pathways <- read_csv('Generic_Models/Recon3D_Rxn_Rules_Pathways.csv')
recon_pathways$SYMBOL <- recon_pathways$Rule
str_detect(unique(recon_pathways$Pathway),"d")
colnames(fva_all)
colnames(recon_pathways)
fva_all <- left_join(fva_all,recon_pathways[,c("Rxn","SYMBOL")])
## Calculate reaction presence for the three glioma models
reaction_presence <- left_join(recon_pathways[,c("Pathway","Rxn")],fva_all[,c("Rxn","Subtype","Flux")])
reaction_presence$Flux[reaction_presence$Flux=="NaN"] <- NA
reaction_presence %>% group_by(Pathway) %>%
  mutate(N_total_rxns = n_distinct(Rxn)) %>%
  group_by(Pathway,Subtype) %>%
  mutate(N_rxns =  n_distinct(Rxn[!is.na(Flux)]),
         Pathway_presence = N_rxns/N_total_rxns) -> reaction_presence
reaction_presence <- unique(reaction_presence[,c("Pathway","N_total_rxns","Subtype","Pathway_presence")])
reaction_presence <- reaction_presence[!is.na(reaction_presence$Subtype),]
#Sort the pathways by the differences between GBM / ((ATS=ODG)/2)
reaction_presence %>% group_by(Pathway)%>%
  mutate(Rank = Pathway_presence[Subtype=="GBM"] - ((Pathway_presence[Subtype=="AST"]+Pathway_presence[Subtype=="ODG"])/2),
        # Rank = abs(Rank)*N_total_rxns
         ) -> reaction_presence
reaction_presence$Pathway_presence <- reaction_presence$Pathway_presence*100
reaction_presence$Pathway_presence_ <- lapply(as.numeric(reaction_presence$Pathway_presence) , function(x) signif.ceiling(x, 3))
reaction_presence$Subtype <- factor(reaction_presence$Subtype,c("GBM","AST","ODG"))
reaction_presence$Subtype_ <-str_c("i",reaction_presence$Subtype)
reaction_presence$Subtype_ <- factor(reaction_presence$Subtype_,c("iGBM","iAST","iODG"))

P_XX_pathways_1 <- ggplot(reaction_presence,aes(x=Subtype_,y=reorder(Pathway,abs(Rank)),
                                     fill=Pathway_presence,label=Pathway_presence_))+
  geom_tile(size=0.4)+#stat = "identity", position=position_dodge())+ #aes(color=Is_Resistant)
  #ggtitle("Drug sensitivity of the anti-melanoma drugs in the Primary PRISM database") +
  geom_text(size=3)+
  scale_fill_gradient2(low = "white",high = "red",guide = "legend",limits = c(0, 100))+
  ylab("") +
  xlab("Subtype models") +theme_classic() +
  theme(axis.text.y = element_text(angle = 0,size = 10.5),
  axis.text.x = element_text(angle = 0,size = 12))+
  guides(fill=guide_legend(title='Reactions per\nPathway (%)')) +
  #theme(legend.position = "left")+
  theme(legend.position=c(-0.75,0.4))
P_XX_pathways_1
max(reaction_presence$N_total_rxns)
max(log10(reaction_presence$N_total_rxns))
reaction_presence_ <- unique(reaction_presence[,c("Pathway","Rank","N_total_rxns")])
P_XX_pathways_2 <- ggplot(reaction_presence_,aes(x=N_total_rxns,y=reorder(Pathway,abs(Rank))))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=N_total_rxns),size=3,nudge_x = 0.3)+
  
  ylab("") +
  xlab("Number of reacions") +theme_classic() +
  scale_x_continuous(trans = "log10")+
  theme(axis.text.y = element_text(angle = 0,size = 0),
        axis.text.x = element_text(angle = 0,size = 12))
P_XX_pathways_2
#P_XX_pathways_final <- cowplot::plot_grid(P_XX_pathways_1, P_XX_pathways_2, ncol = 2,rel_widths = c(6,2.5))
P_XX_pathways_final <- P_XX_pathways_1 + plot_spacer()+ P_XX_pathways_2+ plot_layout(ncol = 3,widths = c(5, -0.4 ,3))#
P_XX_pathways_final
png(filename="Figure/FigS_XX_Pathway_presence.png", units="in", width=10, height=11, res=300)
P_XX_pathways_final
dev.off()
#ggsave('Figure/FigS_XX_Pathwy_presence.png', P_XX_pathways_final,units = 'in',width =11,height = 11,dpi = 300)

# Map drug's deleted reactions to their FVA results
dko_df__  = read_csv('Integrated_Drug_Db/DKO_result_with_Effective_Targets.csv')
sko_df_del <- dko_df__[,c("Drugs","Subtype","DelRxns")]
sko_df_del %>% separate_rows(DelRxns,sep = "; ")-> sko_df_del
# Define subtype-sepfic rxn deletion
sko_df_del %>% group_by(Drugs, DelRxns) %>% mutate(N_subtypes = n_distinct(Subtype)) -> sko_df_del
sko_df_del_spec <- sko_df_del[sko_df_del$N_subtypes<4,]
# focus only on combination drugs
#sko_df_del_spec <- sko_df_del_spec[sko_df_del_spec$Drugs %in% dko_df$Drugs,]
sko_df_del_spec$Subtype <- factor(sko_df_del_spec$Subtype,c("GBM","AST","ODG"))
selected_combinations <- c(#"cannabidiol ;adapalene", "eflornithine ;rifamycin", 
                           #"cannabidiol ;zidovudine",#"fluorouracil ;zidovudine",   
                           #"fluorouracil ;dorzolamide" ,"fluorouracil ;methazolamide", 
                           "fluorouracil ;resveratrol" ,
                          #"fluorouracil ;sulfanilamide"
                           "fluorouracil ;zonisamide" 
                           )
sko_df_del_spec <- sko_df_del_spec[sko_df_del_spec$Drugs %in% selected_combinations,]
PX <- ggplot(sko_df_del_spec,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
               aes(y=DelRxns,x=Subtype))+
  geom_tile(alpha=0.7,stat = 'identity') +  
  theme(axis.text.x = element_text(angle = 90,size = 10),
        axis.text.y = element_text(angle = 0,size = 9))+
  facet_wrap(.~Drugs,scales = "free")
PX
sko_df_del_fva <- left_join(sko_df_del_spec,fva_all,by=c("DelRxns"="Rxn"))
n_distinct(sko_df_del_fva$DelRxns)
sko_df_del_fva$Subtype.y <- factor(sko_df_del_fva$Subtype.y,c("GBM","AST","ODG"))
sko_df_del_fva <- sko_df_del_fva[!str_detect(sko_df_del_fva$Formula,"statin"),]
P_S_combination_rxns <- ggplot(sko_df_del_fva,aes(y=reorder_within(str_c(DelRxns, ": ", Formula),rank,rank_2),#reorder_within(str_c(Rxn, ": ", Formula),rank,rank_2),
                                             alpha=FluxType,x = abs(Flux),
                                       group=Subtype.y,
                                       shape=Reaction,
                                       #color=class,fill=class,
                                       #color=sub_class,fill=sub_class,
                                       size=FluxType))+
  geom_point(position = position_dodge(width = 0.7))+ #
  #geom_linerange(aes(xmin=abs(ymin),xmax=abs(ymax)), size = 0.6, alpha = 0.5,
  #               position = position_dodge(width = 0.7))+
  ggtitle("Uptake flux rate in the three glioma subtype models") +
  #scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'))+
  scale_y_reordered()+
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(2,3))+
  xlab("Absolute flux rate") +
  ylab("Reactions") +theme_classic() +
  facet_grid(Drugs~Subtype.y,scales = "free")+ #Subtype
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 9),
  )
P_S_combination_rxns



## Visualize the fluxes of cholesterol metabolism pathways ("Cholesterol metabolism"             "Squalene and cholesterol synthesis")
fva_all_ <- left_join(fva_all,recon_pathways)
fva_x <- fva_all[str_detect(fva_all_$Formula,"lutamate") & str_detect(fva_all$Formula,"ctadecan"),]
#fva_all <- fva_all[!is.na(fva_all$Rule),]
unique(recon_pathways$Pathway)[str_detect(unique(recon_pathways$Pathway),"holesterol")]
chol_pathways <- c("Cholesterol metabolism","Squalene and cholesterol synthesis")
chol_genes <- unique(recon_pathways_longer$SYMBOL[recon_pathways_longer$Pathway %in% chol_pathways])
fva_chol <- fva_all[fva_all_$Pathway %in% chol_pathways,]
fva_chol <- fva_chol[fva_chol$average_flux>1e-4,]
n_distinct(fva_chol$Rxn[fva_chol$average_flux>1e-4])

P_S_cholesterol <- ggplot(fva_chol,aes(y=reorder_within(str_c(Rxn, ": ", Formula),rank,rank_2),alpha=FluxType,x = abs(Flux),
                       group=Subtype,
                       shape=Reaction,
                       #color=class,fill=class,
                       #color=sub_class,fill=sub_class,
                       size=FluxType))+
  geom_point(position = position_dodge(width = 0.7))+ #
  #geom_linerange(aes(xmin=abs(ymin),xmax=abs(ymax)), size = 0.6, alpha = 0.5,
  #               position = position_dodge(width = 0.7))+
  ggtitle("Uptake flux rate in the three glioma subtype models") +
  #scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'))+
  scale_y_reordered()+
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(2,3))+
  xlab("Absolute flux rate") +
  ylab("Reactions") +theme_classic() +
  facet_grid(.~Subtype,scales = "free")+ #Subtype
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 10),
  )
P_S_cholesterol
ggsave(P_S_cholesterol,filename="Figure/FigS_X_Cholesterol_FVA.png", units="in", width=16, height=12)
P_S_cholesterol <- ggplot(fva_all_[fva_all_$super_class %in% "Lipids and lipid-like molecules",],
                          aes(y=reorder_within(str_c(Rxn, ": ", Formula),rank,rank_2),alpha=FluxType,x = abs(Flux),
                                       group=Subtype,
                                       shape=Reaction,
                                       #color=class,fill=class,
                                       #color=sub_class,fill=sub_class,
                                       size=FluxType))+
  geom_point(position = position_dodge(width = 0.7))+ #
  #geom_linerange(aes(xmin=abs(ymin),xmax=abs(ymax)), size = 0.6, alpha = 0.5,
  #               position = position_dodge(width = 0.7))+
  ggtitle("Uptake flux rate in the three glioma subtype models") +
  #scale_color_manual(values = c('#BB173A','#2F6790','forestgreen'))+
  scale_y_reordered()+
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(2,3))+
  xlab("Absolute flux rate") +
  ylab("Reactions") +theme_classic() +
  facet_grid(.~Subtype,scales = "free")+ #Subtype
  theme(axis.text.x = element_text(angle = 0,size = 11),
        axis.text.y = element_text(angle = 0,size = 10),
  )
P_S_cholesterol
ggsave(P_S_cholesterol,filename="Figure/FigS_X_Cholesterol_FVA.png", units="in", width=16, height=12)

## Rank of cholesterol genes in Depmap
colnames(sko_all)
chol_df <- data.frame(Targets = chol_genes,Drugs="_",Subtype="AST:GBM:ODG",TYPE2="Cholesterol metabolism",Predicted_Essential = "No")

chol_df %>% separate_rows(Subtype,sep=':')%>%distinct()->chol_df
sko_all_chol <- rbind(sko_all,chol_df)

depmap_dep_chol <- depmap_dep_longer_rank_placeholder
drugTypes_2 <- c("Single drug target","Combination drug target","Essential gene", "Anti-brain cancer target","Cholesterol metabolism")

depmap_dep_chol$TYPE2 <- paste0(drugTypes_2,collapse = ":")
depmap_dep_chol %>% separate_rows(TYPE2,sep=':')%>%distinct()->depmap_dep_chol

depmap_dep_chol <- left_join(depmap_dep_chol,sko_all_chol) 


#depmap_dep_chol$TYPE2[depmap_dep_chol$Targets %in% chol_genes] <- "Cholesterol metabolism"
depmap_dep_chol$TYPE2 <- factor(depmap_dep_chol$TYPE2,drugTypes_2)
all_genes <- unique(all_genes)
all_genes_2 <- union(all_genes,chol_genes)

#Keeping only genes in the previous four classes
depmap_dep_chol <- depmap_dep_chol[depmap_dep_chol$Targets %in% all_genes_2,]
depmap_dep_chol <- depmap_dep_chol[depmap_dep_chol$Subtype != "Non-glioma",]
#depmap_dep_chol <- depmap_dep_chol[!(is.na(depmap_dep_chol$Predicted_Essential) & depmap_dep_chol$TYPE2 =="Anti-brain cancer target"), ]
depmap_dep_chol <- depmap_dep_chol[!(is.na(depmap_dep_chol$Predicted_Essential) & depmap_dep_chol$TYPE2 =="Cholesterol metabolism"), ]

#depmap_dep_chol$Predicted_Essential[depmap_dep_chol$TYPE2 =="Cholesterol metabolism"] <- "No"
#depmap_dep_chol$Drugs[depmap_dep_chol$TYPE2 =="Cholesterol metabolism"] <- "_"

#depmap_dep_chol <- depmap_dep_chol[!(is.na(depmap_dep_chol$Predicted_Essential) & depmap_dep_chol$TYPE2 =="Cholesterol metabolism"), ]

depmap_dep_chol %>% group_by(Targets) %>%
  mutate(to_drop = TYPE2 %in% TYPE2[!is.na(Predicted_Essential)]) -> depmap_dep_chol
depmap_dep_chol <- depmap_dep_chol[depmap_dep_chol$to_drop==T,]

depmap_dep_chol <- unique(depmap_dep_chol)
depmap_dep_chol %>% group_by(Targets,Subtype) %>%# 
  mutate(rank=max(rank),median_prob=max(median_prob))%>% distinct() ->depmap_dep_chol
depmap_dep_chol$rank_1 <- 0
depmap_dep_chol <- depmap_dep_chol[!is.na(depmap_dep_chol$median_prob),]
depmap_dep_chol$rank_1[depmap_dep_chol$Subtype=='GBM'] <- depmap_dep_chol$rank[depmap_dep_chol$Subtype=='GBM']

n_distinct(depmap_dep_chol$Targets)
n_distinct(depmap_dep_chol$Drugs)


depmap_drugs<- depmap_dep_chol
#depmap_dep_chol <- distinct(depmap_dep_chol[,!colnames(depmap_dep_chol) %in% c("Drugs","Tested_InVitro")])

depmap_dep_chol <- depmap_dep_chol%>% group_by(Targets,TYPE2) %>% 
  mutate(median_prob_gbm = unique(median_prob[Subtype=='GBM']))

depmap_dep_chol$Predicted_per_subtype <- "No"
depmap_dep_chol$Predicted_per_subtype[!is.na(depmap_dep_chol$Predicted_Essential)] <- "Yes"


shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

# Make font map for the essential genes
depmap_essential_bold <- depmap_dep_chol[depmap_dep_chol$Subtype=='GBM',]  |> 
  select(Targets,rank_1,Predicted_Essential)|>
  distinct() |>
  dplyr::arrange(desc(rank_1))
depmap_essential_bold$Predicted_Essential[depmap_essential_bold$Predicted_Essential=='Yes'] <- "bold"
depmap_essential_bold$Predicted_Essential[depmap_essential_bold$Predicted_Essential=='No'] <- "plain"
depmap_dep_chol$Subtype_n_cells <- factor(depmap_dep_chol$Subtype_n_cells ,c("GBM (n=44)","AST (n=20)","ODG (n=2)"))
variant_genes <- unique(depmap_dep_chol$Targets[depmap_dep_chol$Predicted_per_subtype=="No"])

depmap_dep_chol[,colnames(depmap_dep_chol)!= "TYPE2"] %>% 
  reframe(Targets,median_prob_gbm) %>%
  distinct( )%>%
  dplyr::arrange(-median_prob_gbm)->depmap_dep_longer_genes

n_distinct(depmap_dep_longer_genes$Targets)
genes_sorted_1 = depmap_dep_longer_genes$Targets[1:81]
genes_sorted_2 = depmap_dep_longer_genes$Targets[82:162]
unique(depmap_dep_chol$TYPE2)
unique(depmap_dep_chol$Targets) 

P_XX_depmap_1 <- ggplot(depmap_dep_chol[depmap_dep_chol$Targets %in%  genes_sorted_1,],
                        aes(y=reorder(Targets,median_prob_gbm),
                            fill=TYPE2,
                            alpha=Predicted_per_subtype,
                            x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  #ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the glioma cell lines out of 17386 genes") +
  ggtitle("A") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug targets/ essential genes") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue',"violet" ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue',"violet"))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 9)
  )+
  theme(legend.position = c(0.9,0.35),legend.background = element_blank(),legend.key = element_blank())
P_XX_depmap_1

P_XX_depmap_2 <- ggplot(depmap_dep_chol[depmap_dep_chol$Targets %in%  genes_sorted_2,],
                        aes(y=reorder(Targets,median_prob_gbm),
                            fill=TYPE2,
                            alpha=Predicted_per_subtype,
                            x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=TYPE2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("B") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug targets/ essential genes") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~Subtype_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  scale_fill_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue',"violet" ))+#"#105d46"
  scale_color_manual(values = c("forestgreen",'#ffa100',"#d1b430",'steelblue',"violet"))+
  
  scale_alpha_discrete(range=c(0.3,1))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill='none',color='none',alpha="none")+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        #axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_Essential ),
        axis.text.y = element_text(angle = 0,size = 8.5)
  )
P_XX_depmap_2
P_XX_depmap_final <- cowplot::plot_grid(P_XX_depmap_1, P_XX_depmap_2, ncol = 2,rel_widths = c(5,2.5))
P_XX_depmap_final
ggsave(P_XX_depmap_final,filename="Figure/FigS_XX_DepMap_Rank_Cholesterol.png", units="in", width=13, height=11.8, dpi=300)



## PCA of the selected samples expression data

library(DESeq2)
library(knitr)
library(arrayQualityMetrics)
library(devtools)

#Read TCGA data
fpkm = read_csv('data/TCGA_Rahman2015/TCGA_LGG_GBM.csv')

rowname <- fpkm$`Unnamed: 0`
row.names(fpkm ) <- fpkm$`Unnamed: 0`
fpkm =fpkm[,2:ncol(fpkm)]

#fpkm = as.data.frame(t(fpkm))
colname= gsub('\\-','_',colnames(fpkm) )
colnames(fpkm) <-  substr(colname,1,nchar(colname))

#Read TCGA metadata
meta = read_csv('data/TCGA_TCGBiolinks_metadata_Summary.csv')
meta$barcode = gsub('\\-','_',meta$barcode)
table(meta$WHO_2021)

unk_samples <- meta$barcode[!meta$WHO_2021 %in% c('NA_','ODG_AST')]
#intersect(unk_samples,colname)

fpkm <- fpkm[,colname %in% unk_samples]
#fpkm <- fpkm[,colname %in% meta$barcode]

meta <- meta[meta$barcode %in% unk_samples,]
meta <- meta[meta$barcode %in% colname,]
library("FactoMineR")
library("factoextra")
PC <- PCA(t(fpkm), graph = FALSE)

P_X_PCA_1 <- fviz_pca_ind(PC, axes=c(1,2), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='A PC1 vs PC2')
P_X_PCA_2 <- fviz_pca_ind(PC, axes=c(1,3), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='B) PC1 vs PC3')
P_X_PCA_3 <- fviz_pca_ind(PC, axes=c(2,3), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='C) PC2 vs PC3')

#P_X_PCA_4 <- fviz_pca_ind(PC, axes=c(3,4), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2,
#                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='PC3 vs PC4')
P_X_PCA_Final <- (P_X_PCA_1 ) /( P_X_PCA_2 + P_X_PCA_3 )+
  guides(color=guide_legend(title="WHO 2021\nClassification"))+ plot_layout(guides = "collect")


png(filename="Figure/FigS_X_PCA.png", units="in", width=6, height=7, res=300)
P_X_PCA_Final
dev.off()

# PCA over metabolic expression data 
fpkm_metabolic <- fpkm[rowname%in% recon_pathways_longer$SYMBOL,]
PC_metabolic <- PCA(t(fpkm_metabolic), graph = FALSE)
P_X_PCA_1 <- fviz_pca_ind(PC_metabolic, axes=c(1,2), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='A PC1 vs PC2')
P_X_PCA_2 <- fviz_pca_ind(PC_metabolic, axes=c(1,3), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='B) PC1 vs PC3')
P_X_PCA_3 <- fviz_pca_ind(PC_metabolic, axes=c(2,3), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2, #col.ind="WHO_2021",
                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='C) PC2 vs PC3')

#P_X_PCA_4 <- fviz_pca_ind(PC, axes=c(3,4), col.ind = meta$WHO_2021, pointshape = 20, pointsize = 2,
#                          geom =  c("point"), addEllipses = TRUE, ellipse.type = "confidence",  ellipse.level=0.95, title='PC3 vs PC4')
P_X_PCA_metabolic_Final <- (P_X_PCA_1 ) /( P_X_PCA_2 + P_X_PCA_3 )+
  guides(color=guide_legend(title="WHO 2021\nClassification"))+ plot_layout(guides = "collect")


png(filename="Figure/FigS_X_PCA_metabolic_genes.png", units="in", width=8, height=7, res=300)
P_X_PCA_metabolic_Final
dev.off()

n_distinct(sko_df_genes[sko_df_genes$TYPE2 =="Single Drug Deletion","SYMBOL"])





## Testing using super cript for Figure 2.B
# library(ggtext)
# x <- unique(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',
#                          c("Drugs","Drugs_one_ch")])
# x$Drugs_one_ch_superscript <- str_c(x$Drugs_one_ch,"<sup>2</sup>")
# 
# sko_df_genes_$Drugs_one_ch_superscript <- str_c(sko_df_genes_$Drugs_one_ch,"<sup>2</sup>")
# 
# 
# P2_2_ <- ggplot(sko_df_genes_,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
#                aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank))+
#   geom_tile(alpha=0.7,stat = 'identity',aes(fill=Drugs)) +
#   geom_richtext(aes(label=Drugs_one_ch_superscript,color=Drugs,group=Drugs_one_ch_superscript, size=Predicted_Essential),
#             position="identity" ,label.size = NA)+ #,size=2.5
#   facet_grid(Predicted_Essential~Subtype_,scales = 'free')+#, space='free_y')+
#   
#   force_panelsizes(rows = c(1.2, 2.9),
#                    cols = c(1, 1)) +
#   
#   ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
#   ylab("")+xlab("Number of drugs per gene")+
#   theme_classic() +
#   scale_fill_manual(values=mycols)+
#   scale_color_manual(values=rep("black",33),labels=sort(x$Drugs))+
#   scale_size_manual(values=c(3.2,2))+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 9) ,
#         strip.text = element_markdown(size = 14),
#         legend.text=element_text(size=12),legend.title = element_text(size=13),
#         axis.title = element_text(size = 12))+
#   #axis.text.y= element_markdown(size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0,size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0, 
#   #                           color = ifelse(sko_df_genes_$defensive_industries == "N", "red", "black"))
#   
#   guides(color=guide_legend(title="Drugs",ncol = 1,keywidth = 0.5),
#          fill=guide_legend(override.aes = list(label = sort(x$Drugs_one_ch_superscript))),size="none"
#   ) 
# P2_2_
# P2_final <- cowplot::plot_grid(P2, P2_2_, ncol = 2,rel_widths = c(3.3,5.5))
# P2_final
# ggsave(P2_final,filename="Figure/Figure_2_Essential_Genes_.png", units="in", width=13, height=14,dpi=300)
# 
# ## Testing using stroke for Figure 2.B
# effective_df_keep$Drugs <- firstup(effective_df_keep$Drugs)
# sko_df_genes_indication <- left_join(sko_df_genes_,effective_df_keep[,c("Drugs","Indication")])
# n_distinct(sko_df_genes_indication$Indication)
#                                                                      
# P2_2_ <- ggplot(sko_df_genes_indication,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
#                aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank))+
#   #geom_tile(alpha=0.7,stat = 'identity',aes(fill=Drugs)) +
#   geom_point(alpha=0.7,stat = 'identity',size = 6,aes(color=Drugs,shape=Indication)) +
#   geom_text(aes(label=Drugs_one_ch,group=Drugs_one_ch, size=Predicted_Essential,color=Drugs),position="identity" )+ #,size=2.5
#   facet_grid(Predicted_Essential~Subtype_,scales = 'free')+#, space='free_y')+
#   
#   force_panelsizes(rows = c(1.5, 2.7),
#                    cols = c(1, 1)) +
#   
#   ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
#   ylab("")+xlab("Number of drugs per gene")+
#   theme_classic() +
#   scale_fill_manual(values=mycols)+
#   scale_color_manual(values=rep("black",33),labels=sort(x$Drugs))+
#   scale_size_manual(values=c(4,2.5))+
#   scale_shape_manual(values=c(0,1,2,21,22,23,24,25))+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 9) ,
#         strip.text = element_text(size = 14),
#         legend.text=element_text(size=12),legend.title = element_text(size=13),
#         axis.title = element_text(size = 12))+
#   #axis.text.y= element_markdown(size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0,size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0, 
#   #                           color = ifelse(sko_df_genes_$defensive_industries == "N", "red", "black"))
#   
#   guides(fill=guide_legend(title="Drugs",ncol = 1,keywidth = 0.5),
#          color=guide_legend(override.aes = list(label = sort(x$Drugs_one_ch))),size="none"
#   ) 
# P2_2_
# 
# P2_final <- cowplot::plot_grid(P2, P2_2_, ncol = 2,rel_widths = c(3.3,5.5))
# P2_final
# ggsave(P2_final,filename="Figure/Figure_2_Essential_Genes_.png", units="in", width=13, height=14,dpi=300)
# 
# 
# 
# ### Stroked tile
# P2_2_ <- ggplot(sko_df_genes_indication,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
#                aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank))+
#   geom_tile(alpha=0.7, width=0.95, height=0.95, size=2,stat = 'identity',aes(fill=Drugs,color1=Indication)) +
#   geom_point(alpha=0.8, size=6,stat = 'identity',aes(shape = Drugs),fill="black") +
#   
#   geom_text(aes(label=Drugs_one_ch,group=Drugs,color=Drugs, size=Predicted_Essential,#color=Drugs,
#                 ),position="identity" )+ #,size=2.5,color="black"
#   facet_grid(Predicted_Essential~Subtype_,scales = 'free')+#, space='free_y')+
#   
#   force_panelsizes(rows = c(1.5, 2.7),
#                    cols = c(1, 1)) +
#   
#   ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
#   ylab("")+xlab("Number of drugs per gene")+
#   theme_classic() +
#   scale_fill_manual(values=mycols)+
#   scale_color_manual(aesthetics  ="color1",
#                       values = c("#F8766D", "#FF61CC", "#CD9600", "#7CAE00", "#C77CFF", "#00BE67" ,"#00BFC4", "#00A9FF"))+
#   #scale_color_manual(aesthetics ="color2",values =mycols)+
#   scale_shape_manual(values=seq(1,33),labels=sort(x$Drugs))+
#   #scale_colour_manual(aesthetics = c("color1", "color2"),
#   #                   colours = list(c("#F8766D", "#FF61CC", "#CD9600", "#7CAE00", "#C77CFF", "#00BE67" ,"#00BFC4", "#00A9FF"),
#   #                                  mycols))+
#   #scale_color_manual(values=rep("black",33),labels=sort(x$Drugs))+
#   scale_size_manual(values=c(4,2.5))+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 9) ,
#         strip.text = element_text(size = 14),
#         legend.text=element_text(size=12),legend.title = element_text(size=13),
#         axis.title = element_text(size = 12))+
#   #axis.text.y= element_markdown(size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0,size = 8,face =bold_map$FACE )
#   #axis.text.y = element_text(angle = 0, 
#   #                           color = ifelse(sko_df_genes_$defensive_industries == "N", "red", "black"))
#   
#   guides(fill=guide_legend(title="Drugs",ncol = 1,keywidth = 0.5),
#          shape=guide_legend(override.aes = list(label = sort(x$Drugs_one_ch),title="Drugs",color="black")),
#          size="none"
#   ) 
# P2_2_


library(ggnewscale)
effective_df_keep$Drugs <- firstup(effective_df_keep$Drugs)
sko_df_genes_indication <- left_join(sko_df_genes_,effective_df_keep[,c("Drugs","Indication")])
n_distinct(sko_df_genes_indication$Indication)

sko_df_genes_indication %>% group_by(SYMBOL_,Predicted_Essential,Subtype) %>%mutate(
  SYMBOL_id = str_length(paste0(unique(Drugs),collapse = "; " )),
  Drug_len = str_length(Drugs),n_drugs = n_distinct(Drugs),
  Drug_rank_2 = ntile(str_length(Indication),n_distinct(Drugs)),
) ->sko_df_genes_indication


P2_2 <- ggplot(sko_df_genes_indication,#(sko_df_genes[sko_df_genes$TYPE2=='Single Drug Deletion',],
               aes(y=reorder(SYMBOL_,n_drugs),x=Drug_rank_2))+
  geom_raster(alpha=0.7,stat = 'identity',aes(fill=Drugs)) +

  geom_text(aes(label=Drugs_one_ch,color=Drugs,group=Drugs_one_ch, size=Predicted_Essential),position="identity" )+ #,size=2.5
  facet_grid(Predicted_Essential~Subtype_,scales = 'free')+#, space='free_y')+
  
  force_panelsizes(rows = c(1.5, 2.7),
                   cols = c(1, 1)) +
  
  ggtitle("B") +#  #ggtitle("Predicted drug targets") +# 
  ylab("")+xlab("Number of drugs per gene")+
  theme_classic() +
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=rep("black",33),labels=sort(x$Drugs))+
  scale_size_manual(values=c(4.2,3.3))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 11) ,
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),legend.title = element_text(size=13),
        axis.title = element_text(size = 12))+
  #axis.text.y= element_markdown(size = 8,face =bold_map$FACE )
  #axis.text.y = element_text(angle = 0,size = 8,face =bold_map$FACE )
  #axis.text.y = element_text(angle = 0, 
  #                           color = ifelse(sko_df_genes_$defensive_industries == "N", "red", "black"))
  
  guides(fill=guide_legend(title="Drugs",ncol = 1,keywidth = 0.5),
         #color=guide_legend(override.aes = list(label = sort(x$Drugs_one_ch))),
         size="none"
  ) +
  # start a new scale
  new_scale_color() +
  geom_tile(alpha=0.7, width=0.9, height=0.9, size=1.25,stat = 'identity',aes(color=Indication),fill="NA") +
  scale_color_manual(values = c("#F8766D", "#00A9FF", "red", "#00BE67", "#CD9600", "#7CAE00", "#C77CFF" ,"#00BFC4"))
  
  #geom_rect(size=1, fill=NA,#labels= unique(sko_df_genes_indication$Indication),
  #           aes(colour=Indication,xmin=Drug_rank - 0.5, xmax=Drug_rank + 0.5,
  #                ymin=reorder(SYMBOL_,n_drugs) , ymax=reorder(SYMBOL_,n_drugs) ))
  #scale_color_identity(legend = FALSE)

P2_final <- cowplot::plot_grid(P2, P2_2, ncol = 2,rel_widths = c(3.3,5.2))
P2_final
ggsave(P2_final,filename="Figure/Figure_2_Essential_Genes_2.png", units="in", width=12.5, height=12.5,dpi=300)

### Model deleted reactions and rea
sko_delrxn <- read_csv("./Integrated_Drug_Db/SKO_Drugs_Del_Rxns.csv")
