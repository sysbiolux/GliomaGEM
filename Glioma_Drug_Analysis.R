## Glioma drug analysis 
# 06/03/2022

#1- Number of clinical trials by phase and completion ( https://beta.clinicaltrials.gov/)
#2- Main targets (Drug Repurposing Hub)
#3- Targets in approved or on-going trials in melanoma (Drug Repurposing Hub, ReDO_Trials_DB)
#5- Rank of the target genes in CRISPR data on melanoma cell lines (metastasis /normal /unknown) (DepMap 22Q1)
#6- Rank of the drugs using viability assay on melanoma cell lines (metastasis /normal /unknown) (primary PRISM database 19Q4)


## What are the databases used:

## ReDO_Trials_DB 
# "is a curated database of active clinical trials investigating the use of non-cancer drugs as potential cancer treatments"
# https://www.anticancerfund.org/en/redo-trials-db

## Drug Repurposing Hub
# A well-curated databases of approved drugs and their targets.

## DepMap 22Q1
# A pan-cancer crispr screen over xx cell lines and xx genes

## Primary PRISM database
# A high-throuput cell vibaility screen on xx cancer cell lines and xx drugs

# clinical_targets = the targets of approved anti-melanoma drugs + on-going clinical trials in melanoma
# Is any of the predicted drugs' targets are in clinical_targets

setwd("./GliomGEM")
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
library(ggrepel)
library(DT)
library(vroom)
library(readxl)
library(data.table)
library(tidyverse)
library(dplyr)
library(clipr)
library(devtools)
library(gridExtra)
library(grid)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(knitr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(rrcov)
library(vroom)
library(tibble)
library(lemon)
library(ggplotify)
library(patchwork)
library(readr)
library(tidytext)
library(vegalite)
library(lamisc)
library(scales)

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#3- Targets in approved or on-going trials in melanoma (Drug Repurposing Hub, ReDO_Trials_DB)

approved_drugs <- drh_longer[str_detect(drh_longer$indication,'melanoma'),]
approved_drugs <- as.data.frame(approved_drugs[!is.na(approved_drugs$Drugs),])
approved_drugs[is.na(approved_drugs$target),] # Check approved drugs with missing targets
approved_drugs[approved_drugs$Drugs=='fotemustine','target'] <- 'TXNRD1' # fotemustine was missing the target, filled from moa
antimelanoma <- approved_drugs[,c("Drugs","target" ,"Approved_Anticancer")]
antimelanoma_2 <- read.csv('./Anti-Melanoma_Drugs.txt')
setdiff(antimelanoma_2$Aldesleukin,antimelanoma$Drugs)

colnames(antimelanoma) <- c("Drugs","Targets" ,"Tested_InVitro")
anti_NOS <- drh_longer[str_detect(drh_longer$moa,'nitric oxide synthase inhibitor'),c("Drugs","target")]
anti_NOS <- anti_NOS[!is.na(anti_NOS$Drugs),]
colnames(anti_NOS) <- c("Drugs","Targets" )
anti_NOS$MOA <- 'Nitric oxide synthase inhibitor'
anti_NOS$Source <- 'Drug Repurpusing Hub'
anti_NOS_2 <- read_csv('./Anti-NO_Drugs_PMC7912608.csv')
anti_NOS_2$Targets <- ''

#anti_NOS_ <- rbind(distinct(anti_NOS[,c('Drugs','MOA','Source')]),anti_NOS_2[,c('Drugs','MOA','Source')])
#write_clip(anti_NOS_)

anti_NOS <- rbind(anti_NOS,anti_NOS_2[,colnames(anti_NOS)])
anti_NOS$Drugs <- str_replace(anti_NOS$Drugs,'\\-',' ')
anti_NOS$Tested_InVitro <- "NO-based drug"
anti_NOS$Drugs <- str_to_lower(anti_NOS$Drugs)

# Remove manually a mislabelled enzyme as drug
anti_NOS <- anti_NOS[anti_NOS$Drugs !='aminomethyltransferase', ] 

intersect(approved_drugs$target,drugs_df_longer$Targets) # None of the predicted drugs share targets with approved anti-melanoma drugs
intersect(approved_drugs$target,drugs_df_longer_all$All_Targets) # None of the predicted drugs share targets with approved anti-melanoma drugs

## 
ReDO_Trials_DB <- read_delim('ReDO_Trials_DB.txt',delim='\t')
ReDO_Trials_DB$Drug_INN <- str_to_lower(ReDO_Trials_DB$Drug_INN)
intersect(ReDO_Trials_DB$Drug_INN,drugs_df_longer$Drugs) # None of the predicted drugs share targets with approved anti-melanoma drugs

# On-going clinical trials for cancer using the predicted drugs
drugs_ongoing_trials <- left_join(drugs_df_drh, ReDO_Trials_DB,by=c('Drugs'='Drug_INN'))
drugs_ongoing_trials <- drugs_ongoing_trials[!is.na(drugs_ongoing_trials$`NCT Number`),]
colnames(drugs_ongoing_trials)
drugs_ongoing_trials <- drugs_ongoing_trials[,c("Drugs","Tested_InVitro",'Approved_Anticancer' ,"NCT Number","Phase","Status","Enrollment" ,"Interventions" , "Conditions", "Cancer_Group", "Cancer_Type",
                                                "Stage","Controlled","Multi-Arm", "Primary-EP")]
drugs_ongoing_trials$Enrollment<- as.numeric(drugs_ongoing_trials$Enrollment)
# Correcting the phase of a study https://ascopubs.org/doi/10.1200/JCO.2020.38.15_suppl.e22081
drugs_ongoing_trials$Phase[drugs_ongoing_trials$`NCT Number`=='NCT03972748'] <- 'Phase 2'
drugs_ongoing_trials$Phase[drugs_ongoing_trials$Phase=='Not available/Missing'] <- 'Missing'
#drugs_ongoing_trials[drugs_ongoing_trials$Phase=='Missing',] <- NA
drugs_ongoing_trials <- drugs_ongoing_trials[!drugs_ongoing_trials$Phase=='Missing',]
drugs_ongoing_trials$Phase <- factor(drugs_ongoing_trials$Phase,c("Phase 4","Phase 3","Phase 2/3" ,"Phase 2","Phase 1/2", "Phase 1"))
write_csv(drugs_ongoing_trials,'Figure/Non-cancer_trials.csv')
p4 <- ggplot(drugs_ongoing_trials,aes(y=Drugs,x=Enrollment,color=Cancer_Group,shape=Status,
                                      size=Tested_InVitro
))+
  geom_point()+
  scale_shape_manual(values = seq(1,10,1)) +
  scale_x_continuous(trans='log10')+
  theme_classic()+
  xlab("Number of participants") +
  ylab("Non-cancer drugs") +
  #guides(alpha=guide_legend(title="Finished/\nStill recruiting?"))+
  #ylab("Total number of predicted essential genes") +
  facet_grid(Phase~.,scales = 'free_y')+
  scale_size_discrete(range=c(3,7))+
  ggtitle("Repurposed non-cancer drugs in cancer clinical trials")+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))
p4

ggsave('Figure/Fig4_Clinical_trials_in_cancer_for_noncancer_drugs.png', p4,units = 'in',
       width = 9,height = 9,dpi = 300)


## Trials of NOS based drugs
# On-going clinical trials for cancer using the predicted drugs
drugs_ongoing_trials <- left_join(anti_NOS, ReDO_Trials_DB,by=c('Drugs'='Drug_INN'))
drugs_ongoing_trials <- drugs_ongoing_trials[!is.na(drugs_ongoing_trials$`NCT Number`),]
colnames(drugs_ongoing_trials)
drugs_ongoing_trials <- drugs_ongoing_trials[,c("Drugs","Tested_InVitro" ,"NCT Number","Phase","Status","Enrollment" ,"Interventions" , "Conditions", "Cancer_Group", "Cancer_Type",
                                                "Stage","Controlled","Multi-Arm", "Primary-EP")]
drugs_ongoing_trials <- distinct(drugs_ongoing_trials)
drugs_ongoing_trials$Enrollment<- as.numeric(drugs_ongoing_trials$Enrollment)
# Correcting the phase of a study https://ascopubs.org/doi/10.1200/JCO.2020.38.15_suppl.e22081
drugs_ongoing_trials$Phase[drugs_ongoing_trials$`NCT Number`=='NCT03972748'] <- 'Phase 2'
drugs_ongoing_trials$Phase[drugs_ongoing_trials$Phase=='Not available/Missing'] <- 'Missing'
#drugs_ongoing_trials[drugs_ongoing_trials$Phase=='Missing',] <- NA
drugs_ongoing_trials <- drugs_ongoing_trials[!drugs_ongoing_trials$Phase=='Missing',]
drugs_ongoing_trials$Phase <- factor(drugs_ongoing_trials$Phase,c("Phase 4","Phase 3","Phase 2/3" ,"Phase 2","Phase 1/2", "Phase 1"))
write_csv(drugs_ongoing_trials,'Figure/Non-cancer_trials.csv')
p4_ <- ggplot(drugs_ongoing_trials,aes(y=Drugs,x=Enrollment,color=Cancer_Group,shape=Status,
                                       size=Tested_InVitro))+
  geom_point()+
  scale_shape_manual(values = seq(1,10,1)) +
  scale_x_continuous(trans='log10')+
  theme_classic()+
  xlab("Number of participants") +
  ylab("Non-cancer drugs") +
  #guides(alpha=guide_legend(title="Finished/\nStill recruiting?"))+
  #ylab("Total number of predicted essential genes") +
  facet_grid(Phase~.,scales = 'free_y')+
  scale_size_discrete(range=c(3,7))+
  ggtitle("Repurposed non-cancer drugs in cancer clinical trials")+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))
p4_

ggsave('Figure/Fig4__Clinical_trials_NOS_drugs.png', p4,units = 'in',
       width = 9,height = 9,dpi = 300)



## Rank in primary PRISM database
# Define glioma cell lines and make a separate column for it
depmap_cellinfo = read.csv("./DepMap_22Q1/sample_info.csv")

depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID','cell_line_name',"Subtype",
                                      'primary_or_metastasis','sex','depmap_public_comments','CCLE_Name')]
depmap_cellinfo <- depmap_cellinfo[str_detect(depmap_cellinfo$Subtype,'Melanoma'),]
depmap_cellinfo$primary_or_metastasis[depmap_cellinfo$primary_or_metastasis==""] <- "Uncategorized "
depmap_cellinfo$sex[depmap_cellinfo$sex==""] <- "Uncategorized "
depmap_cellinfo$primary_or_metastasis <- factor(depmap_cellinfo$primary_or_metastasis,c("Primary",'Metastasis',"Uncategorized "))
prism_pr <- read_csv('./PRISM_Repurposing_19Q4/primary-screen-replicate-collapsed-logfold-change.csv')
colnames(prism_pr)[1] <- 'DepMap_ID'
prism_pr <- prism_pr[prism_pr$DepMap_ID %in% depmap_cellinfo$DepMap_ID,]
prism_pr <- as.data.frame(prism_pr)
row.names(prism_pr) <- prism_pr$DepMap_ID

prism_pr_drugs <- read_csv('./PRISM_Repurposing_19Q4/primary-screen-replicate-treatment-info.csv')
prism_pr_drugs <- distinct(prism_pr_drugs[,c('broad_id','name')])
prism_pr_drugs

## Define the rank of the sko drugs in the glioma cell lines (overall rank not by subtype)
prism_pr %>% pivot_longer(cols = !DepMap_ID,names_to = "broad_conc",values_to = "Log2FoldChange") -> prism_pr_longer
prism_pr_longer %>% separate(broad_conc,sep = '::',into = c('broad_id','conc')) -> prism_pr_longer
prism_pr_longer <- left_join(prism_pr_longer,prism_pr_drugs)

# selecting only melanoma cell lines
prism_pr_longer <- prism_pr_longer[prism_pr_longer$DepMap_ID%in% depmap_cellinfo$DepMap_ID,]

N_celllines = n_distinct(prism_pr_longer$DepMap_ID)
N_drugs <- n_distinct(prism_pr_longer$name)
n_distinct(prism_pr_longer$broad_id)
n_distinct(prism_pr_longer$name)

prism_pr_longer <- prism_pr_longer[,colnames(prism_pr_longer)!='broad_id']
prism_pr_longer <- prism_pr_longer[!is.na(prism_pr_longer$Log2FoldChange),]

prism_pr_longer <- left_join(prism_pr_longer,depmap_cellinfo)

prism_pr_longer$name <- str_replace(prism_pr_longer$name,'\\-',' ')
prism_pr_longer$name[prism_pr_longer$name =='icatibant acetate'] <- 'icatibant'
prism_pr_longer$name[prism_pr_longer$name =='LGX818'] <- 'encorafenib'
prism_pr_longer$name[prism_pr_longer$name =='MEK162'] <- 'binimetinib'

prism_pr_longer$conc <- lapply(as.numeric(prism_pr_longer$conc) , function(x) signif.ceiling(x, 3))
prism_pr_longer$name_conc <- str_c(prism_pr_longer$name, " (", prism_pr_longer$conc,")")
n_distinct(prism_pr_longer$name_conc)

prism_pr_longer %>% group_by(name_conc,primary_or_metastasis)%>% # Subtype
  mutate(median_logfc = median(Log2FoldChange,na.rm = T),
         FoldChange= (1-(2^Log2FoldChange))*100,
         median_FoldChange = median(FoldChange,na.rm = T))%>%
  group_by(primary_or_metastasis) %>%
  mutate( n_cells= n_distinct(DepMap_ID), # number of cell lines by subtype
          primary_n_cells =  str_c(primary_or_metastasis," (n=", n_cells,")"))  %>%
  group_by(DepMap_ID)%>%
  dplyr::arrange(median_FoldChange) -> prism_pr_longer_rank

prism_pr_longer_rank %>% group_by(name_conc,sex)%>% # Subtype
  mutate(median_FoldChange_sex = median(FoldChange,na.rm = T))%>%
  group_by(sex) %>%
  mutate( n_cells_sex= n_distinct(DepMap_ID), # number of cell lines by sex
          sex_n_cells =  str_c(sex," (n=", n_cells_sex,")"))  -> prism_pr_longer_rank

#rank_in_group(.)-k
#prism_pr_longer_rank$rank <- ntile(prism_pr_longer_rank$median_FoldChange,
#                                            n_distinct(prism_pr_longer_rank$name))
prism_pr_longer_rank %>% group_by(primary_or_metastasis) %>% 
  mutate(rank =  ntile(-median_FoldChange,n_distinct(name_conc))) -> prism_pr_longer_rank
prism_pr_longer_rank$name <- str_to_lower(prism_pr_longer_rank$name)

setdiff(drugs_df$Drugs,prism_pr_longer_rank$name)
setdiff(antimelanoma$Drugs,prism_pr_longer_rank$name)

n_distinct(prism_pr_longer_rank$median_FoldChange)/3

max(prism_pr_longer_rank$rank)
unique(prism_pr_longer_rank$rank)

antimelanoma <- approved_drugs[,c("Drugs","target" ,"Approved_Anticancer","Is_metabolic")]
colnames(antimelanoma) <- c("Drugs","Targets" ,"Tested_InVitro","Is_metabolic")

setdiff(str_to_lower(anti_NOS$Drugs),prism_pr_longer_rank$name)

unique(prism_pr_longer_rank$name)
prism_pr_longer_rank$name[str_detect(prism_pr_longer_rank$name,"S Isopropylisothiourea")]


anti_NOS$Is_metabolic <- "NA"
anti_NOS$Is_predicted <- "NA"
antimelanoma$Is_predicted <- "NA"
x <- drugs_df_longer[,c("Drugs","All_Targets","Tested_InVitro","Is_metabolic",'Is_predicted' )]
colnames(x)[2] <- 'Targets'
drugs_antimelanoma <- rbind(x,
                            antimelanoma,anti_NOS[,c("Drugs","Targets","Tested_InVitro","Is_metabolic",'Is_predicted' )])

drugs_antimelanoma$Tested_InVitro[drugs_antimelanoma$Tested_InVitro=='Yes'] <- 'Approved anti-melanoma'
drugs_antimelanoma$Drugs <- str_to_lower(drugs_antimelanoma$Drugs)
colnames(drugs_df)
prism_pr_longer_rank <- left_join(prism_pr_longer_rank,
                                  distinct(drugs_antimelanoma[,c("Drugs" ,"Tested_InVitro")]),by=c('name'='Drugs')) 

## Define resistant cell lines to approved anti-melanoma drugs
prism_pr_longer_rank$name_conc <- firstup(prism_pr_longer_rank$name_conc)
prism_resist <- prism_pr_longer_rank[prism_pr_longer_rank$Tested_InVitro=='Approved anti-melanoma',]

prism_resist <- prism_resist[!is.na(prism_resist$Tested_InVitro),]
prism_resist %>% group_by(DepMap_ID) %>% # rank cell lines per resistance
  mutate(median_cell = mean(FoldChange,na.rm = T)) -> prism_resist
min_value <- min(prism_resist$FoldChange)
max_value <- max(prism_resist$FoldChange)
prism_resist$FoldChange_ <- lapply(as.numeric(prism_resist$FoldChange) , function(x) signif.ceiling(x, 3))
prism_resist$Is_Resistant <- 'No'
prism_resist$Is_Resistant[prism_resist$FoldChange< -50] <- 'Yes'
prism_resist$Is_Resistant <- factor(prism_resist$Is_Resistant,c('Yes','No'))
p5_resist <- ggplot(prism_resist,aes(x=name_conc,y=reorder(cell_line_name,desc(median_cell)),
                                     fill=FoldChange,label=FoldChange_))+
  geom_tile(aes(color=Is_Resistant),size=0.4)+#stat = "identity", position=position_dodge())+
  ggtitle("Drug sensitivity of the anti-melanoma drugs in the Primary PRISM database") +
  geom_text(size=4)+
  scale_fill_gradientn(colours = c("cyan", "white", "red"), #na.value = na.value.forplot,
                       values = scales::rescale(c(max_value, 0, min_value)))+
  ylab("Cell lines") +
  xlab("Drugs with conc in µM") +theme_classic() +
  scale_color_manual(values=c('black','white'))+
  #facet_grid(Type~.,scales = 'free')+ # ,space="free"
  theme(axis.text.x = element_text(angle = 90,size = 12),
        axis.text.y = element_text(angle = 0,size = 12),
        panel.background = element_rect(fill = "grey92", colour = NA))+
  guides(fill=guide_legend(title='Reduction in viability\nrelative to DMSO (%)'))
p5_resist
ggsave('Figure/SupplementaryFigureS11_New_resistance.png', p5_resist,units = 'in',width = 8,height = 10,dpi = 300)
ggsave('Figure/SupplementaryFigureS11_New_resistance.pdf', p5_resist,units = 'in',width = 8,height = 10,dpi = 300)


####
prism_pr_longer_rank %>% group_by(DepMap_ID) %>% # rank cell lines per resistance in anti-melanoma
  filter(!is.na(Tested_InVitro)) %>% 
  mutate(Antimelanoma_Resistance = mean(FoldChange[Tested_InVitro=='Approved anti-melanoma']
                                        ,na.rm = T)) -> prism_resist_db

prism_resist_db <- prism_resist_db[,c('name_conc','Tested_InVitro','FoldChange','cell_line_name','CCLE_Name'
                                      ,'Antimelanoma_Resistance',
                                      'primary_or_metastasis','sex')]
celllines_of_interest <- c('A101D_SKIN','A2058_SKIN','A375_SKIN','C32_SKIN','CJM_SKIN','COLO679_SKIN','COLO741_SKIN','COLO783_SKIN','COLO792_SKIN','COLO800_SKIN','COLO829_SKIN','G361_SKIN','HMCB_SKIN','HS294T_SKIN','HS695T_SKIN','HS852T_SKIN','HS936T_SKIN','HS939T_SKIN','HS944T_SKIN','HT144_SKIN','IGR1_SKIN','IGR37_SKIN','IGR39_SKIN','IPC298_SKIN','K029AX_SKIN','LOXIMVI_SKIN','MALME3M_SKIN','MDAMB435S_SKIN','MELHO_SKIN','MELJUSO_SKIN','MEWO_SKIN','RPMI7951_SKIN','RVH421_SKIN','SH4_SKIN','SKMEL1_SKIN','SKMEL24_SKIN','SKMEL28_SKIN','SKMEL30_SKIN','SKMEL31_SKIN','SKMEL3_SKIN','SKMEL5_SKIN','UACC257_SKIN','UACC62_SKIN','WM115_SKIN','WM1799_SKIN','WM2664_SKIN','WM793_SKIN','WM88_SKIN','WM983B_SKIN')
#prism_resist_db <- prism_resist_db[prism_resist_db$CCLE_Name %in%celllines_of_interest, ]
write_csv(prism_resist_db,'./Figure/PRISM_resistance_data.csv')

# top_resistant cell lines by selecting cellines with more than 50% proliferation
prism_resist_cells <- prism_resist_db[prism_resist_db$FoldChange <=-50,]
prism_resist_cells <- prism_resist_cells[prism_resist_cells$Tested_InVitro=='Approved anti-melanoma',]
n_distinct(prism_resist_cells$cell_line_name)
# Concatentae the restant drugs tegether
prism_resist_cells %>% group_by(cell_line_name) %>% 
  mutate(Drug= str_split(name_conc," \\(",simplify=T )[,1]) %>%
  group_by(cell_line_name) %>%
  summarize(Resistance = paste0(sort(firstup(Drug)),collapse =  ' & ')) ->prism_resist_cells
prism_resist_cells$Resistance
prism_resist_cells$Resistance <- str_replace(prism_resist_cells$Resistance ,
                                             "Binimetinib & Cobimetinib & Dacarbazine & Trametinib",
                                             "Binimetinib & Cobimetinib &\nDacarbazine & Trametinib")
prism_resist_db <- left_join(prism_resist_db, prism_resist_cells)
prism_resist_db<- prism_resist_db %>% group_by(name_conc) %>%
  mutate(score_resistant_cells = median(FoldChange[cell_line_name %in% prism_resist_cells$cell_line_name]))
# Fold change for the selected points (resistant cell lines)
prism_resist_db$FoldChange_sel <- NA
prism_resist_db$FoldChange_sel[!is.na(prism_resist_db$Resistance)] <-prism_resist_db$FoldChange[!is.na(prism_resist_db$Resistance)]
prism_resist_db$cell_line_name_sel <- NA
prism_resist_db$cell_line_name_sel[!is.na(prism_resist_db$Resistance)] <-prism_resist_db$cell_line_name[!is.na(prism_resist_db$Resistance)]

#Create a custom color scale
library(RColorBrewer)
prism_resist_db$Tested_InVitro[prism_resist_db$Tested_InVitro %in% c("Selected","Not selected")] <- "Predicted drug"

prism_resist_color <- prism_resist_db  |> 
  select(name_conc,score_resistant_cells,Tested_InVitro)|>
  distinct() |>
  dplyr::arrange(desc(-score_resistant_cells))
#myColors <- brewer.pal(3,"Dark2")
#names(myColors) <- unique(prism_resist_color$Tested_InVitro)
#myColors <- as.data.frame(myColors)
#myColors$Tested_InVitro <- rownames(myColors)
myColors <- data.frame(myColors=c("forestgreen", "#D95F02" ,"#7570B3"),
                       Tested_InVitro=c("Predicted drug", "Approved anti-melanoma","NO-based drug"))

prism_resist_color <- left_join(prism_resist_color,myColors)
prism_resist_db$Tested_InVitro <- factor(prism_resist_db$Tested_InVitro,
                                         c("Predicted drug","Approved anti-melanoma","NO-based drug"))
p5_resist_ <- ggplot(prism_resist_db,aes(y=reorder(name_conc,desc(-score_resistant_cells)),#p5_resist_ <- 
                                         #,
))+
  geom_violin(aes(x=FoldChange))+
  geom_point(aes(shape=cell_line_name_sel,x=FoldChange_sel,color=Resistance))+#,color=Resistance #
  geom_point(aes(x=score_resistant_cells),color='black',shape=11)+#,color=Resistance #
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.5)+
  ggtitle("Sensitivity analysis of the anti-melanoma resistant\ncell lines in the Primary PRISM database") +
  xlab("Reduction in viability relative to DMSO (%)") +
  ylab("Drugs with conc in µM") +theme_classic() +
  scale_shape_manual(values=seq(1,10,1))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_markdown(angle = 0,size = 12, color =prism_resist_color$myColors ),
        panel.background = element_rect( colour = NA))+ #fill = "grey92"
  guides(color=guide_legend(title='Resistance'),
         shape=guide_legend(title='Resistant cell lines'))+
  geom_dotplot(aes(x=FoldChange,fill = Tested_InVitro), alpha = 0,binwidth = 9, key_glyph = draw_key_text) +#
  scale_fill_manual(values = c('Predicted drug'="red",
                               "Approved anti-melanoma"="#D95F02" ,
                               "NO-based drug"="#7570B3"),guide=guide_legend(title = "Drug type",
                                                                             override.aes = list(alpha = 1, size = 8,color=c("forestgreen","#D95F02","#7570B3"))))+
  theme(legend.position = c(0.15, 0.35),legend.background = element_blank(),legend.key = element_blank())
#guides(fill = ) #+
#p5_resist_
ggsave('Figure/SupplementaryFigureS11_New_resistance_AllDrugs.png', p5_resist_,units = 'in',width = 10,height = 12,dpi = 300)
ggsave('Figure/SupplementaryFigureS11_New_resistance_AllDrugs.pdf', p5_resist_,units = 'in',width = 10,height = 12,dpi = 300)

colnames(prism_resist_db)

#################

#prism_pr_longer <- left_join(prism_pr_longer,depmap_cellinfo[,c('DepMap_ID','Subtype')])


length(unique(prism_pr_longer$DepMap_ID))

# Calculate the n percentile of the median rank
prism_pr_longer_rank$PRISM_Primary <- ntile(prism_pr_longer_rank$median_FoldChange,100)
Ranking_df_1 <- distinct(prism_pr_longer_rank[!is.na(prism_pr_longer_rank$Tested_InVitro),c('name','PRISM_Primary')])

length(unique(prism_pr_longer_rank$name))
prism_pr_longer_rank %>% group_by(name_conc,primary_or_metastasis) %>%# 
  mutate(rank=max(rank))%>% distinct() ->prism_pr_longer_rank
colnames(prism_pr_longer_rank)
prism_pr_longer_rank <- distinct(prism_pr_longer_rank[,c('name','median_FoldChange','median_FoldChange_sex',
                                                         'rank',
                                                         "Tested_InVitro",'conc','primary_or_metastasis',"sex",
                                                         'name_conc','primary_n_cells','sex_n_cells')])
#prism_pr_longer_rank$clinical_phase[is.na(prism_pr_longer_rank$clinical_phase)] <- 'Uncategorized '
prism_pr_longer_rank <- na.omit(prism_pr_longer_rank)
max(prism_pr_longer_rank$rank)
# Sort by primary cell lines rank

prism_pr_longer_rank$rank_1 <- 0
prism_pr_longer_rank$rank_1[prism_pr_longer_rank$primary_or_metastasis=='Metastasis'] <- prism_pr_longer_rank$rank[prism_pr_longer_rank$primary_or_metastasis=='Metastasis']
prism_pr_longer_rank$rank_1[prism_pr_longer_rank$primary_or_metastasis!='Metastasis'] <- 5000

#prism_pr_longer_rank$Tested_InVitro <- factor(prism_pr_longer_rank$Tested_InVitro,
#                                              c("Selected","Not selected","Approved Anti-melanoma","NO-based drug"))
prism_pr_longer_rank$Tested_InVitro[prism_pr_longer_rank$Tested_InVitro =='Selected'] <- "Predicted and tested"
prism_pr_longer_rank$Tested_InVitro[prism_pr_longer_rank$Tested_InVitro =='Not selected'] <- "Predicted and untested"

#prism_pr_longer_rank$Tested_InVitro <- factor(prism_pr_longer_rank$Tested_InVitro,
#                                              c("Predicted and tested","Predicted and untested","Approved Anti-melanoma","NO-based drug"))

prism_pr_longer_rank$DrugType <- prism_pr_longer_rank$Tested_InVitro
prism_pr_longer_rank$DrugType[prism_pr_longer_rank$DrugType %in% c("Predicted and tested","Predicted and untested")] <- "Predicted drug"
prism_pr_longer_rank$DrugType <- factor(prism_pr_longer_rank$DrugType,
                                        c("Predicted drug","Approved anti-melanoma","NO-based drug"))
prism_pr_longer_rank$name_conc <- firstup(prism_pr_longer_rank$name_conc)
p5 <- ggplot(prism_pr_longer_rank,aes(y=reorder(name_conc,desc(rank_1)),
                                      fill=DrugType,
                                      x = median_FoldChange))+
  geom_bar(stat = "identity", position=position_dodge())+ # ,aes(alpha=Tested_InVitro)
  geom_text(aes(label=rank),size=3.7)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("Rank of the predicted drugs in the Primary PRISM\ndatabase in the melanoma cell lines out of 4604 drugs") +
  xlab("Median reduction in viability relative to DMSO (%)") +
  ylab("Drugs and tested concentration (µM)") +theme_classic() +
  facet_grid(.~primary_n_cells)+
  scale_fill_brewer(palette="Dark2")+
  guides(fill=guide_legend(title="Drug class"))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))
p5
#ggsave('Figure/Fig5_Primary_PRISM_database_rank.png', p5,units = 'in',width = 11,height = 10,dpi = 300)
ggsave('Figure/SupplementaryFigureS11_New.png', p5,units = 'in',width = 11,height = 12,dpi = 300)
ggsave('Figure/SupplementaryFigureS11_New.pdf', p5,units = 'in',width = 11,height = 12,dpi = 300)

# ## By sex
## No differences in the cell lines using the sex status
# p5_sex <- ggplot(prism_pr_longer_rank,aes(y=reorder(name_conc,desc(rank_1)),
#                                       fill=DrugType,
#                                       x = (1-median_FoldChange_sex)*100))+
#   geom_bar(stat = "identity", position=position_dodge())+ # ,aes(alpha=Tested_InVitro)
#   #geom_text(aes(label=rank),size=3.7)+
#   geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
#   ggtitle("Rank of the predicted drugs in the Primary PRISM\ndatabase in the melanoma cell lines out of 4604 drugs") +
#   xlab("Median reduction in viability relative to DMSO (%)") +
#   ylab("Drugs and tested concentration (µM)") +theme_classic() +
#   facet_grid(.~sex_n_cells)+
#   scale_fill_brewer(palette="Dark2")+
#   guides(fill=guide_legend(title="Drug class"))+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12))
# p5_sex
# #ggsave('Figure/SupplementaryFigureS11_New_by_sex.png', p5_sex,units = 'in',width = 11,height = 10,dpi = 300)

setdiff(ordered_drugs,prism_pr_longer_rank$name)
table(distinct(prism_pr_longer_rank[,c('name','DrugType')])['DrugType'])
prism_pr_longer_rank %>% group_by(name_conc) %>% 
  mutate(median_FoldChange_meta =median(median_FoldChange[primary_or_metastasis=='Metastasis'] )) ->prism_pr_longer_rank

prism_pr_longer_rank_ <-prism_pr_longer_rank[prism_pr_longer_rank$median_FoldChange_meta >= 50 ,]
p5_ <- ggplot(prism_pr_longer_rank_,aes(y=reorder(name_conc,desc(rank_1)),
                                        fill=DrugType,
                                        x = median_FoldChange))+
  geom_bar(stat = "identity", position=position_dodge())+ # ,aes(alpha=Tested_InVitro)
  geom_text(aes(label=rank,color=DrugType),size=3.7,nudge_x = 10)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("B")+
  #ggtitle("Rank of the predicted drugs in the Primary PRISM\ndatabase in the melanoma cell lines out of 4604 drugs") +
  xlab("Median reduction in viability relative to DMSO (%)") +
  ylab("Drugs and tested concentration (µM)") +theme_classic() +
  facet_grid(.~primary_n_cells)+
  #scale_alpha_discrete(range=c(1,0.5))+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  guides(fill=guide_legend(title="Drug class"),color='none')+
  theme(plot.title= element_text(angle = 0,size = 18,face="bold"),
        axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_text(angle = 0,size = 12))
p5_


###

depmap_dep = read_csv('./DepMap_22Q1/CRISPR_gene_dependency.csv')
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
unique(depmap_dep_longer$depmap_public_comments)
depmap_dep_longer$Resistance <- str_replace(depmap_dep_longer$depmap_public_comments,"Drug resistance: ","")

depmap_dep_longer %>% 
  group_by(primary_or_metastasis) %>%
  mutate( n_cells= n_distinct(DepMap_ID), # number of cell lines by subtype
          primary_n_cells =  str_c(primary_or_metastasis," (n=", n_cells,")"))  %>%
  group_by(Targets) %>%# 
  mutate(median_prob_per_gene=median(Dependancy_Probability,na.rm = T))%>%
  group_by(Targets,primary_or_metastasis,primary_n_cells,median_prob_per_gene) %>%# 
  summarise(median_prob=median(Dependancy_Probability,na.rm = T))%>%
  dplyr::arrange(median_prob)->depmap_dep_longer_rank

depmap_dep_longer_rank %>% group_by(primary_n_cells) %>% 
  mutate(rank =  ntile(-median_prob,n_distinct(Targets))) -> depmap_dep_longer_rank

essential_genes <- c("PTPMT1","PGS1","CRLS1","SGMS1","SPTLC1","SPTLC2","SPTLC3","KDSR","CMPK1","TXNRD1","CAD","DHODH","UMPS","GUK1","FDFT1","SQLE","LSS","CYP51A1","MSMO1","EBP","TM7SF2","DHCR7","HMGCR","MVK","MVD","PMVK","NSDHL","ACACA","LCAT","LIPA","FASN","HSD17B4","ANPEP","SLC27A4","SLC7A5")
essential_df <- data.frame(Targets= essential_genes, Predicted_essential="Yes")
NOS_genes <- c("NOS1","NOS2","NOS3","ARG1","ARG2","ARG3","ASL","ASS1","PDE1A","PDE4A","PDE5A",
               "GUCY1A2", "GUCY1A3", "GUCY1B2" , "GUCY1B3")
NOS_df <- data.frame(Targets= NOS_genes, Predicted_essential="No")
essential_df <- rbind(essential_df,NOS_df)
intersect(essential_genes,antimelanoma$Targets)

# Also add non-metabolic targets
#drugs_antimelanoma$Is_Metabolic <- 'Metabolic'
#drugs_antimelanoma_nonmet <- drugs_df_longer[drugs_df_longer$Targets!=drugs_df_longer$Non_metabolic_targets,
#                                             c("Drugs","Non_metabolic_targets" ,"Tested_InVitro")]
colnames(drugs_antimelanoma) 
#colnames(drugs_antimelanoma_nonmet) <- c("Drugs" , "Targets","Tested_InVitro")

drugs_antimelanoma_2 <- merge(x = drugs_antimelanoma[drugs_antimelanoma$Tested_InVitro!='NO-based drug',],
                              y = essential_df, by = "Targets", all = TRUE)
#drugs_antimelanoma_2 <- merge(x = drugs_antimelanoma_2, y = NOS_df, by = "Targets", all = TRUE)

drugs_antimelanoma_2$Predicted_essential[is.na(drugs_antimelanoma_2$Predicted_essential)] <- 'No'
#depmap_dep_longer_rank$Target_type <- "Drug target"

drugs_antimelanoma_2$Gene_type <- 'Drug target'
drugs_antimelanoma_2$Gene_type[drugs_antimelanoma_2$Tested_InVitro %in% "Approved anti-melanoma"] <-"Anti-melanoma target"
drugs_antimelanoma_2$Gene_type[drugs_antimelanoma_2$Predicted_essential %in% "Yes" &
                                 is.na(drugs_antimelanoma_2$Tested_InVitro) ] <-"Non-druggable essential"
drugs_antimelanoma_2$Gene_type[drugs_antimelanoma_2$Targets %in% NOS_genes ] <- "NO-related genes"
drugs_antimelanoma_2 <- unique(drugs_antimelanoma_2)

drugs_antimelanoma_2$Is_metabolic[drugs_antimelanoma_2$Is_metabolic=='Yes'] <- "Metabolic"
drugs_antimelanoma_2$Is_metabolic[drugs_antimelanoma_2$Is_metabolic=='No'] <- "Non-metabolic"
drugs_antimelanoma_2$Is_predicted[is.na(drugs_antimelanoma_2$Is_predicted)] <- "NA"
drugs_antimelanoma_2$Gene_type_2 <- drugs_antimelanoma_2$Gene_type


idx <- drugs_antimelanoma_2$Gene_type %in% c("Anti-melanoma target",'Drug target')
drugs_antimelanoma_2$Gene_type_2[idx] <- str_c(
  drugs_antimelanoma_2$Is_metabolic[idx]," ", drugs_antimelanoma_2$Gene_type[idx] )
# Define main targets 
# All anti-melanoma drugs are main target

# drugs_antimelanoma_2 <- left_join(drugs_antimelanoma_2,drh_longer[,c('Drugs','target','Approved_Anticancer')],
#                                   by=c("Drugs"="Drugs",'Targets'='target'))
# drugs_antimelanoma_2$Gene_type[drugs_antimelanoma_2$Gene_type=='Predicted drug target' &
#                                  !is.na(drugs_antimelanoma_2$Approved_Anticancer)] <- "Predicted drug main target"
# drugs_antimelanoma_2$Gene_type[drugs_antimelanoma_2$Gene_type=='Predicted drug target' &
#                                    is.na(drugs_antimelanoma_2$Approved_Anticancer)] <- "Predicted drug off-target"

#drugs_antimelanoma_2$Target_type <- NA
#drugs_antimelanoma_2$Target_type[!is.na(drugs_antimelanoma_2$Approved_Anticancer)] <- 'Main target'
#drugs_antimelanoma_2$Target_type[drugs_antimelanoma_2$Gene_type=='Predicted drug target' &
#                                   is.na(drugs_antimelanoma_2$Approved_Anticancer)] <- 'Off-target'

depmap_dep_longer_rank <- left_join(drugs_antimelanoma_2  ,depmap_dep_longer_rank)
unique(depmap_dep_longer_rank$Gene_type_2)

#depmap_dep_longer_rank$Gene_type <- factor(depmap_dep_longer_rank$Gene_type,
#                                         c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
#                                           "Approved anti-melanoma target","NO-related genes"))
#depmap_dep_longer_rank$Gene_type <- factor(depmap_dep_longer_rank$Gene_type,
#                                         c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
#                                           "Approved anti-melanoma target","NO-related genes"))

## Approved_or_Predicted
# Approved 
# Predicted

## Prediction
# Antimelanoma target
# Drug target
# Essential gene


depmap_dep_longer_rank %>% group_by(Targets,primary_or_metastasis) %>%# 
  mutate(rank=max(rank))%>% distinct() ->depmap_dep_longer_rank
depmap_dep_longer_rank$rank_1 <- 0
depmap_dep_longer_rank <- depmap_dep_longer_rank[!is.na(depmap_dep_longer_rank$median_prob),]
depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$primary_or_metastasis=='Metastasis'] <- depmap_dep_longer_rank$rank[depmap_dep_longer_rank$primary_or_metastasis=='Metastasis']
depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$primary_or_metastasis!='Metastasis'] <- 18000
#depmap_dep_longer_rank$Tested_InVitro[depmap_dep_longer_rank$Tested_InVitro %in% c("Selected","Not selected")] <- "Predicted drug target"
n_distinct(depmap_dep_longer_rank$Targets)
n_distinct(depmap_dep_longer_rank$Drugs)


#depmap_dep_longer_rank$Gene_type <- factor(depmap_dep_longer_rank$Gene_type,
#                                           c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
#                                             "Approved anti-melanoma target","NO-related genes"))

depmap_drugs<- depmap_dep_longer_rank
depmap_dep_longer_rank_ <- distinct(depmap_dep_longer_rank[,!colnames(depmap_dep_longer_rank) %in% c("Drugs","Tested_InVitro")])

depmap_dep_longer_rank_ <- depmap_dep_longer_rank_%>% group_by(Targets) %>% 
  mutate(median_prob_meta = unique(median_prob[primary_or_metastasis=='Metastasis']))



shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

# Make font map for the essential genes
depmap_essential_bold <- depmap_dep_longer_rank_[depmap_dep_longer_rank_$primary_or_metastasis=='Metastasis',]  |> 
  select(Targets,rank_1,Predicted_essential)|>
  distinct() |>
  dplyr::arrange(desc(rank_1))
depmap_essential_bold$Predicted_essential[depmap_essential_bold$Predicted_essential=='Yes'] <- "bold"
depmap_essential_bold$Predicted_essential[depmap_essential_bold$Predicted_essential=='No'] <- "plain"

p6 <- ggplot(depmap_dep_longer_rank_[,],aes(y=reorder(Targets,desc(rank_1)),alpha=Is_predicted,
                                            fill=Gene_type_2,
                                            #alpha=Target_type,
                                            x = (median_prob*100)))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=Gene_type_2),size=3,nudge_x = 9)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the melanoma cell lines out of 17386 genes") +
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug Targets") +theme_classic() +
  scale_x_continuous(limits=c(0,110))+
  facet_grid(.~primary_n_cells)+
  #scale_fill_brewer(palette="Dark2")+
  
  #scale_fill_manual(values = c('#105d46',"#22AF84","#d1b430", "#D95F02" ,"#7570B3"))+#"#105d46"
  #scale_color_manual(values = c("#105d46","#22AF84","#d1b430", "#D95F02" ,"#7570B3"))+
  
  #scale_alpha_discrete(range=c(1,0.3))+
  #scale_x_continuous(limits=c(-10,110))+
  #scale_x_continuous(trans = shift_trans(-3))+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  scale_alpha_manual(values = c(1,0.5,1))+
  theme(axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_markdown(angle = 0,size = 9, face =depmap_essential_bold$Predicted_essential ),
        #axis.text.y = element_text(angle = 0,size = 9)
  )+
  theme(legend.position = c(0.55,0.6),legend.background = element_blank(),legend.key = element_blank())
#p6
#ggsave('Figure/Fig6_DepMap_rank.png', p6,units = 'in',width = 12,height = 14,dpi = 300)
ggsave('Figure_2/SupplementaryFigureS10_New.png', p6,units = 'in',width = 10,height = 18,dpi = 300)
ggsave('Figure_2/SupplementaryFigureS10_New.pdf', p6,units = 'in',width = 10,height = 18,dpi = 300)

depmap_dep_longer_rank_ <- depmap_dep_longer_rank_[depmap_dep_longer_rank_$median_prob_meta>=0.5,]
depmap_essential_bold_ <- depmap_essential_bold[depmap_essential_bold$Targets %in%depmap_dep_longer_rank_$Targets, ]
p6_ <- ggplot(depmap_dep_longer_rank_,aes(y=reorder(Targets,desc(rank_1)),
                                          fill=Gene_type_2,
                                          x = median_prob*100))+
  geom_bar(stat = "identity", position=position_dodge())+
  geom_text(aes(label=rank,color=Gene_type_2),size=3.7,nudge_x = 12)+
  geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
  #ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the melanoma cell lines out of 17386 genes") +
  ggtitle("A")+
  xlab("Median of the dependency probability per gene (%)") +
  ylab("Drug Targets") +theme_classic() +
  scale_x_continuous(limits=c(0,115))+
  facet_grid(.~primary_n_cells)+
  guides(fill=guide_legend(title="Gene class"),color='none')+
  #scale_fill_manual(values = c("#105d46","#ffcc33", "#D95F02"))+
  #scale_color_manual(values = c("#105d46","#ffcc33", "#D95F02"))+
  #scale_x_continuous(trans = shift_trans(-3))+
  theme(plot.title= element_text(angle = 0,size = 18,face="bold"),
        axis.text.x = element_text(angle = 0,size = 12),
        axis.text.y = element_markdown(angle = 0,size = 11,face =depmap_essential_bold_$Predicted_essential ),
  )
p6_/ p5_
ggsave('Figure_depmap for main targets//Melanoma_paper_Figure5New.pdf', p6_/ p5_,units = 'in',width = 11,height = 10,dpi = 300)
ggsave('Figure_depmap for main targets//Melanoma_paper_Figure5New.png', p6_/ p5_,units = 'in',width = 11,height = 9,dpi = 300)

### Rank drugs per median gene probability
depmap_dep_longer_rank %>% group_by(Drugs) %>% 
  mutate(median_prob_per_drug= median(median_prob_per_gene)) ->depmap_dep_longer_rank
depmap_dep_longer_rank$Drugs_2 <- depmap_dep_longer_rank$Drugs
depmap_dep_longer_rank$Drugs_2[depmap_dep_longer_rank$Gene_type_2=='NO-related genes'] <- 'NO-related genes'
depmap_dep_longer_rank$Drugs_2[depmap_dep_longer_rank$Gene_type_2=='Non-druggable essential'] <- 'Non-druggable essential'

n_distinct(depmap_dep_longer_rank$Drugs)
depmap_dep_longer_rank_drug <- unique(depmap_dep_longer_rank[,c('Drugs_2','Targets','median_prob_per_gene','median_prob_per_drug',
                                                                'Is_metabolic','Is_predicted',"Gene_type")])
sorted_drugs <- unique(depmap_dep_longer_rank_drug$Drugs_2[order(depmap_dep_longer_rank_drug$median_prob_per_drug,decreasing = TRUE)])
depmap_dep_longer_rank_drug$Drugs_2 <- factor(depmap_dep_longer_rank_drug$Drugs_2 ,sorted_drugs)
colors <- rep("brown1",length(sorted_drugs))
colors[sorted_drugs %in% depmap_dep_longer_rank$Drugs[depmap_dep_longer_rank$Tested_InVitro=='Approved anti-melanoma']] <- "cyan"
  colors[sorted_drugs %in% depmap_dep_longer_rank_drug$Drugs_2[depmap_dep_longer_rank_drug$Drugs_2=='NO-related genes']] <- "magenta"
    colors[sorted_drugs %in% depmap_dep_longer_rank_drug$Drugs_2[depmap_dep_longer_rank_drug$Drugs_2=='Non-druggable essential']] <- "darkorange"
      
    write.xlsx(depmap_dep_longer_rank_drug,"Depenedency_per_drug.xlsx")
    p6_by_drugs <- ggplot(depmap_dep_longer_rank_drug,
                          aes(y=reorder(Targets,desc(-median_prob_per_gene)),alpha=Is_predicted,
                              fill=Is_metabolic,
                              #alpha=Target_type,
                              x = (median_prob_per_gene*100)))+
      geom_bar(stat = "identity", position=position_dodge())+
      #geom_text(aes(label=rank,color=Is_metabolic),size=3,nudge_x = 9)+
      geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
      ggtitle("DepMap dependency probability for the predicted drug targets\nand essential genes in the melanoma cell lines out of 17386 genes") +
      xlab("Median of the dependency probability per gene (%)") +
      ylab("Drug Targets") +theme_classic() +
      scale_x_continuous(limits=c(0,110))+
      facet_wrap(.~Drugs_2,ncol = 6,scales = "free_y")+
      guides(fill=guide_legend(title="Gene class"),color='none')+
      scale_alpha_manual(values = c(1,0.5,1))+
      theme(axis.text.x = element_text(angle = 0,size = 12),
            axis.text.y = element_markdown(angle = 0,size = 9)#, face =depmap_essential_bold$Predicted_essential ),
            #axis.text.y = element_text(angle = 0,size = 9)
      )+theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_markdown(fill = "colors"))
    #theme(legend.position = c(0.55,0.6),legend.background = element_blank(),legend.key = element_blank())
    #p6
    #ggsave('Figure/Fig6_DepMap_rank.png', p6,units = 'in',width = 12,height = 14,dpi = 300)
    ggsave('Figure_depmap for main targets//SupplementaryFigureS10_New_by_drugs.png', p6_by_drugs,units = 'in',width = 14,height = 14,dpi = 300)
    ggsave('Figure_depmap for main targets//SupplementaryFigureS10_New_by_drugs.pdf', p6_by_drugs,units = 'in',width = 14,height = 14,dpi = 300)
    
    
    ## Map target genes to pathways
    depmap_rank <- left_join(depmap_dep_longer_rank,recon_summ)
    depmap_rank$Pathway[is.na(depmap_rank$Pathway)] <- "Non-metabolic pathways"
    
    depmap_rank %>% group_by(Pathway) %>% mutate(rank_pathway =mean(rank_1)) ->depmap_rank
    depmap_rank <- distinct(depmap_rank[,c("Pathway","Targets","rank_pathway","Predicted_essential","median_prob",
                                           "Gene_type","primary_or_metastasis",'primary_n_cells')])
    depmap_rank[depmap_rank$primary_or_metastasis!="Metastasis","Targets"] <- ""
    
    p7 <- ggplot(depmap_rank,aes(y=reorder(Pathway,desc(rank_pathway)),
                                 color=Gene_type,label=Targets,shape=Predicted_essential,
                                 x = median_prob*100))+
      geom_point(alpha=0.9)+#position = position_dodge2(width = 1,reverse = F)
      geom_text_repel(size=3)+
      #scale_fill_brewer(palette="Dark2")+
      geom_vline(xintercept=50, linetype="dashed", color = "green", size=1)+
      ggtitle("Pathway analysis of the DepMap dependency probability for\nthe drug targets and essential genes in the melanoma cell lines") +
      xlab("Median of the dependency probability per gene (%)") +
      ylab("Pathways of the drug targets and essential genes") +theme_classic() +
      #scale_x_continuous(limits=c(-10,110))+
      scale_x_continuous(trans = shift_trans(-3))+
      scale_color_manual(values = c("#105d46","#2bdba6","#ffcc33", "#D95F02" ,"#7570B3"))+
      
      facet_grid(.~primary_n_cells)+
      theme(axis.text.x = element_text(angle = 0,size = 12),
            axis.text.y = element_text(angle = 0,size = 10))
    p7
    #ggsave('Figure/Fig7_DepMap_Pathways.png', p7,units = 'in',width = 12,height = 10,dpi = 300)
    ggsave('Figure/SupplementaryFigureS12_New.png', p7,units = 'in',width = 12,height = 10,dpi = 300)
    ggsave('Figure/SupplementaryFigureS12_New.pdf', p7,units = 'in',width = 12,height = 10,dpi = 300)
    
    
    
    ## Plot the gene depeneendanyc by cdrugs
    library(ggplot2)
    
    depmap_drugs <- depmap_drugs |> group_by(Drugs) |> mutate(sum_prob_drug = sum(median_prob))
    depmap_drugs <- depmap_drugs[!is.na(depmap_drugs$Drugs),]
    
    depmap_drugs <- depmap_drugs[order(depmap_drugs$median_prob,decreasing = T),]
    depmap_drugs$Drugs <- firstup(depmap_drugs$Drugs)
    
    p8 <- ggplot(depmap_drugs,aes(#y=#reorder_within(Drugs,sum_prob_drug,median_prob),
      y=reorder(Drugs,sum_prob_drug),
      fill=Gene_type,alpha=median_prob,
      x = median_prob*100))+
      geom_bar(stat='identity',color="white",size=0.2) +
      #scale_y_reordered()+
      ggtitle("Rank of the predicted drugs based on the DepMap dependency\nprobability of the targets in the melanoma cell lines") +
      xlab("Sum of the median dependency probability (%) of the drug targets") +
      ylab("Drugs") +theme_classic() +
      facet_grid(.~primary_n_cells)+
      theme(axis.text.x = element_text(angle = 0,size = 12),
            axis.text.y = element_text(angle = 0,size = 10))+
      scale_fill_manual(values = c("#105d46","#9e1b42", "#D95F02"))+
      guides(alpha=guide_legend(title="Target median\ndependency probability (%)"))+
      scale_alpha_continuous(range = c(0.3,1))
    ggsave('Figure/Fig8_Drug_Rank_by_Target_Dependancy.png', p8,units = 'in',width = 12,height = 10,dpi = 300)
    ggsave('Figure/Fig8_Drug_Rank_by_Target_Dependancy.pdf', p8,units = 'in',width = 12,height = 10,dpi = 300)
    
    
    ## Drug resistant cell line analysis in DepMap
    depmap_dep_longer$Resistance[depmap_dep_longer$Resistance==''] <- 'Wildtype'
    depmap_dep_longer$cell_line <- str_split(depmap_dep_longer$CCLE_Name,'_',simplify = T)[,1]
    resistant_cells <- unique(depmap_dep_longer$CCLE_Name[depmap_dep_longer$Resistance!="Wildtype"])
    wt_cells <- unique(str_split(resistant_cells,'_',simplify = T)[,1])
    depmap_resist <- depmap_dep_longer[depmap_dep_longer$cell_line %in% wt_cells,]
    depmap_resist$CellType <- "Resistant"
    depmap_resist$CellType[str_count(depmap_resist$CCLE_Name,'_')==1] <- 'Wildtype'
    unique(depmap_resist$cell_line_name)
    
    # Rank genes in each cell line and take the relative difference
    depmap_resist %>% 
      group_by(cell_line,Targets) %>%
      mutate(median_resist=median(Dependancy_Probability[CellType=='Resistant'],na.rm = T),
             relative_resistance =median_resist- Dependancy_Probability[CellType=='Wildtype'])%>%
      group_by(Targets) %>%
      mutate(Rank_1 =median(Dependancy_Probability[CellType=='Resistant'],na.rm = T))%>%
      dplyr::arrange(Rank_1)->depmap_resist
    
    depmap_resist %>% group_by(cell_line) %>% 
      mutate(rank =  ntile(-relative_resistance,n_distinct(Targets))) -> depmap_resist
    
    depmap_resist <- left_join(drugs_antimelanoma_2  ,depmap_resist)
    
    depmap_resist %>% group_by(Targets,cell_line) %>%# 
      mutate(rank=max(rank))%>% distinct() ->depmap_resist
    
    #depmap_dep_longer_rank$rank_1 <- 0#
    #depmap_dep_longer_rank <- depmap_dep_longer_rank[!is.na(depmap_dep_longer_rank$median_prob),]
    #depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$primary_or_metastasis=='Metastasis'] <- depmap_dep_longer_rank$rank[depmap_dep_longer_rank$primary_or_metastasis=='Metastasis']
    #depmap_dep_longer_rank$rank_1[depmap_dep_longer_rank$primary_or_metastasis!='Metastasis'] <- 18000
    #depmap_dep_longer_rank$Tested_InVitro[depmap_dep_longer_rank$Tested_InVitro %in% c("Selected","Not selected")] <- "Predicted drug target"
    n_distinct(depmap_dep_longer_rank$Targets)
    n_distinct(depmap_dep_longer_rank$Drugs)
    
    
    #depmap_resist$Gene_type <- factor(depmap_resist$Gene_type,
    #                                           c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
    #                                             "Approved anti-melanoma target","NO-related genes"))
    
    depmap_drugs<- depmap_resist
    depmap_resist <- distinct(depmap_resist[,!colnames(depmap_resist) %in% c("Drugs","Tested_InVitro")])
    
    depmap_resist <- depmap_resist[!is.na(depmap_resist$Dependancy_Probability),]
    
    shift_trans = function(d = 0) {
      scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
    }
    #depmap_resist$Gene_type <- unfactor(depmap_resist$Gene_type)
    
    # Make font map for the essential genes
    depmap_essential_bold <- depmap_resist[,c("Targets","Rank_1","Predicted_essential")]|>
      distinct() |>
      dplyr::arrange(Rank_1)
    depmap_essential_bold <- depmap_essential_bold[order(depmap_essential_bold$Rank_1,decreasing = F),]
    
    depmap_essential_bold$Predicted_essential[depmap_essential_bold$Predicted_essential=='Yes'] <- "bold"
    depmap_essential_bold$Predicted_essential[depmap_essential_bold$Predicted_essential=='No'] <- "plain"
    
    #Create a custom color scale
    library(RColorBrewer)
    library(varhandle)
    depmap_resist$Rank_1 <- as.numeric(depmap_resist$Rank_1)
    levels(depmap_resist$Gene_type)
    depmap_resist_color <- depmap_resist[,c("Targets","Gene_type","Rank_1")]|>
      distinct() |>
      dplyr::arrange(Rank_1)
    
    depmap_resist_color <- depmap_resist_color[order(depmap_resist_color$Rank_1,decreasing = F),]
    #myColors <- brewer.pal(n_distinct(depmap_resist_color$Gene_type),"Dark2")
    #names(myColors) <- unique(depmap_resist_color$Gene_type)
    #myColors <- as.data.frame(myColors)
    #myColors$Gene_type <- rownames(myColors)105d46
    myColors <- data.frame(myColors=c("#21de04","#2bdba6","#ffcc33", "#ed9b0e" ,"#cf68ed"),
                           Gene_type=c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
                                       "Approved anti-melanoma target","NO-related genes"))
    #levels(depmap_resist_color$Gene_type))
    depmap_resist_color <- left_join(depmap_resist_color,myColors)
    depmap_resist_color <- depmap_resist_color[! (depmap_resist_color$Targets%in%c('DHODH','SQLE') &
                                                    depmap_resist_color$Gene_type=="Predicted drug off-target"),]
    n_distinct(depmap_resist_color$Targets)
    n_distinct(depmap_resist$Targets)
    table(depmap_resist_color$Targets)
    depmap_resist$cell_line_name <- str_replace(depmap_resist$cell_line_name,'-','')
    unique(depmap_resist$Resistance)
    depmap_resist$Resistance <- factor(depmap_resist$Resistance ,
                                       c("Wildtype","Dabrafenib" ,"Dabrafenib and Trametinib",
                                         "Dabrafenib and Roxadustat", "SCH772984"),)
    depmap_resist <- depmap_resist[order(depmap_resist$Dependancy_Probability,decreasing = T),]
    depmap_resist$Gene_type <- factor(depmap_resist$Gene_type,c("Predicted drug main target","Predicted drug off-target","Predicted non-druggable essential",
                                                                "Approved anti-melanoma target","NO-related genes"))
    p6_resist <- ggplot(depmap_resist,aes(y=reorder(Targets,Rank_1),
                                          #reorder_within(Targets,Rank_1,Dependancy_Probability),
                                          x = (Dependancy_Probability*100)))+
      geom_bar(stat='identity',aes(fill = Resistance), alpha = 1,position = position_identity() ) +
      geom_point(aes(shape=cell_line_name),size=2)+ #color=Resistance,
      #scale_alpha_discrete(range=c(0.4,1))+
      #geom_text_repel(aes(label = Gene_type, segment.color = Gene_type)) +
      scale_y_reordered()+
      geom_vline(xintercept=50, linetype="dashed", color = "green", size=0.5)+
      ggtitle("Sensitivity analysis using DepMap dependency probability for the predicted drug targets and\nessential genes in the drug-resistant melanoma cell lines") +
      xlab("Dependency probability per gene (%)") +
      ylab("Drug Targets") +#theme_classic() +
      scale_x_continuous(limits=c(0,110))+
      #geom_boxplot(aes(fill = Gene_type), alpha = 0, key_glyph = draw_key_text) +
      #guides(fill = guide_legend(override.aes = list(alpha = 0, size = 8))) +
      facet_grid(.~cell_line)+
      scale_shape_manual(values=seq(1,12,1))+
      guides(fill = guide_legend(reverse = T))+
      theme(axis.text.x = element_text(angle = 0,size = 12),
            axis.text.y = element_markdown(angle = 0,size = 9, fill =depmap_resist_color$myColors,
                                           #box.color = 'black',
                                           face =depmap_essential_bold$Predicted_essential 
            )  ) +   
      geom_point(aes(color = Gene_type), alpha = 0, key_glyph = draw_key_text) +
      scale_color_manual(values = c("#21de04","#2bdba6","#ffcc33", "#ed9b0e" ,"#cf68ed"))+
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 8))) #+
    #theme(axis.text.y = element_text(color = scales::hue_pal()(5), face = 2))
    
    #p6_resist
    
    ggsave('Figure/SupplementaryFigureS10_New_resistance.pdf', p6_resist,units = 'in',width = 12,height = 14,dpi = 300)
    ggsave('Figure/SupplementaryFigureS10_New_resistance.png', p6_resist,units = 'in',width = 12,height = 14,dpi = 300)
    
    ## Compare the viability from PharmacoDB
    pharmacodb_predicted <- read_csv('./PharmacoDB/PharmacoDB_Predicted_drugs.csv')
    pharmacodb_approved <- read_csv('./PharmacoDB/PharmacoDB_Anti_Melanoma_drugs.csv')
    pharmacodb_nos <- read_csv('./PharmacoDB/PharmacoDB_NO_based_drugs.csv')
    pharmacodb_predicted$DrugType <- 'Predicted drug'
    pharmacodb_approved$DrugType <- 'Approved anti-melanoma'
    pharmacodb_nos$DrugType <- 'NO-based drug'
    
    pharmacodb <- rbind(pharmacodb_predicted,pharmacodb_approved,pharmacodb_nos)
    n_distinct(pharmacodb$`drug name`)
    
    # Seelct melanoma cell lines
    pharmacodb_cells <- read_csv('PharmacoDB/cell_annotation_table-2022-12-09T10_24_25-05_00.csv')
    intersect(pharmacodb$`cell name`,depmap_cellinfo$cell_line_name)
    pharmacodb <- pharmacodb[pharmacodb$`cell name`%in% depmap_cellinfo$cell_line_name,]
    n_distinct(pharmacodb$`drug name`)
    unique(pharmacodb$`drug name`)
    
    pharmacodb %>% group_by(DrugType,`drug name`) %>% 
      summarise(Median_IC50 = median(IC50,na.rm = T),
                Std_IC50 = sd(IC50,na.rm = T)) ->pharmacodb_summ
    
    unique(pharmacodb[pharmacodb$`dataset name`%in%c('GDSC1000','GDSC2000'),'drug name'])
    unique(pharmacodb[pharmacodb$`dataset name`%in%c('CTRPv2'),'drug name'])
    
    
    
    ### Integrating PRISM secondary, GDSC1000, GDSC2000, CTRPv2 for IC50 only
    ## Drug, cell line, IC50
    
    ### PRISM secondary ##
    prism_sc <- read_csv('./Drug_Response_Databases/secondary-screen-dose-response-curve-parameters.csv')
    colnames(prism_sc)[2] <- 'DepMap_ID'
    prism_sc <- prism_sc[prism_sc$passed_str_profiling==TRUE,]
    prism_sc <- as.data.frame(prism_sc)
    
    prism_sc <- prism_sc[,c('name','ic50','DepMap_ID','ccle_name')]
    #
    depmap_cellinfo = read.csv("./DepMap_22Q1/sample_info.csv")
    depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID','cell_line_name')] # ,"Subtype",'primary_or_metastasis','sex','depmap_public_comments','CCLE_Name'
    #depmap_cellinfo <- depmap_cellinfo[str_detect(depmap_cellinfo$Subtype,'Melanoma'),]
    prism_sc <- left_join(prism_sc,depmap_cellinfo)
    prism_sc <- prism_sc[,c('name','cell_line_name','ic50')]
    
    ### GDSC1000
    gdsc1000 <- read_excel('./Drug_Response_Databases/GDSC1_fitted_dose_response_24Jul22.xlsx')
    gdsc1000$IC50 <- exp(gdsc1000$LN_IC50)
    colnames(gdsc1000)
    gdsc1000<- gdsc1000[,c("DRUG_NAME","CELL_LINE_NAME","IC50")]
    
    ### GDSC2000
    gdsc2000 <- read_excel('./Drug_Response_Databases/GDSC2_fitted_dose_response_24Jul22.xlsx')
    gdsc2000$IC50 <- exp(gdsc2000$LN_IC50)
    colnames(gdsc2000)
    gdsc2000<- gdsc2000[,c("DRUG_NAME","CELL_LINE_NAME","IC50")]
    
    # # CTRPv2
    # ctrp <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt',delim='\t')
    # colnames(ctrp)
    # ctrp <- ctrp[,c("master_cpd_id","experiment_id" ,"apparent_ec50_umol")] # area_under_curve
    # ctrp_cpd <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt',delim='\t')
    # ctrp_cll <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt',delim='\t')
    # ctrp_exp <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt',delim='\t')
    # ctrp_cpd <- ctrp_cpd[,c('master_cpd_id','cpd_name')]
    # ctrp_cll <- ctrp_cll[,c('master_ccl_id','ccl_name')]
    # ctrp_exp <- ctrp_exp[,c('experiment_id','master_ccl_id')]
    # ctrp_cll <- left_join(ctrp_cll,ctrp_exp) # Joining, by = "master_ccl_id"
    # 
    # ctrp <- left_join(ctrp,ctrp_cpd) # Joining, by = "master_cpd_id"
    # ctrp <- left_join(ctrp,ctrp_cll) # Joining, by = "experiment_id"
    # ctrp <- ctrp[,c('cpd_name','ccl_name','apparent_ec50_umol')]
    
    ### Genentech Cell Line Screening Initiative
    gcsi <- read_delim('./Drug_Response_Databases/data_5_Genentech_Cell_Line_Screening_Initiative_(gCSI).csv',delim = ',')
    gcsi <- gcsi[,c('Perturbagen','Cell_Line','IC50')]
    
    ### Merge the 4 databases
    colnames(prism_sc) <- c('Drugs','cell_line_name','IC50')
    colnames(gdsc1000) <- c('Drugs','cell_line_name','IC50')
    colnames(gdsc2000) <- c('Drugs','cell_line_name','IC50')
    colnames(gcsi) <- c('Drugs','cell_line_name','IC50')
    
    #colnames(ctrp) <- c('Drugs','cell_line_name','IC50')
    
    prism_sc$Database <-'PRISM Secondary'
    gdsc1000$Database <-'GDSC1000'
    gdsc2000$Database <-'GDSC2000'
    gcsi$Database <-'gCSI'
    
    
    
    Merged_ic50 <- rbind(prism_sc,gdsc1000,gdsc2000,gcsi)
    Merged_ic50 <- distinct(Merged_ic50)
    Merged_ic50 <- Merged_ic50[!is.na(Merged_ic50$IC50),]
    Merged_ic50 <- Merged_ic50[Merged_ic50$IC50!='Inf',]
    
    Merged_ic50$Drugs[Merged_ic50$Drugs =='LGX818'] <- 'encorafenib'
    Merged_ic50$Drugs[Merged_ic50$Drugs =='MEK162'] <- 'binimetinib'
    
    Merged_ic50_drugs <-unique(Merged_ic50[,c('Database','Drugs')])
    
    write_csv(Merged_ic50,'./Merged_IC50_databases.csv')
    
    # 
    antimelanoma_2$Tested_InVitro <- 'Approved anti-melanoma'
    drugs_all <- rbind(drugs_df_longer[,c("Drugs" ,"Tested_InVitro")],antimelanoma[,c("Drugs" ,"Tested_InVitro")],
                       antimelanoma_2,anti_NOS[,c("Drugs" ,"Tested_InVitro")])
    drugs_all$Tested_InVitro[drugs_all$Tested_InVitro=='Yes'] <- 'Approved anti-melanoma'
    drugs_all$Tested_InVitro[drugs_all$Tested_InVitro %in% c('Selected', "Not selected")]  <- 'Predicted drug'
    drugs_all <- distinct(drugs_all)
    
    drugs_all$Drugs <- str_to_lower(drugs_all$Drugs)
    
    Merged_ic50$Drugs <- str_to_lower(Merged_ic50$Drugs )
    n_distinct(Merged_ic50$Drugs)
    n_distinct(Merged_ic50$cell_line_name)
    
    #Merged_ic50$Drugs <- str_replace(Merged_ic50$Drugs,'-',' ')
    
    drugs_ic50 <- left_join(drugs_all,Merged_ic50) # select our drugs in the mergedd databse of IC50
    
    # select melanoma cell lines 
    depmap_cellinfo = read.csv("./DepMap_22Q1/sample_info.csv")
    
    depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID','cell_line_name',"Subtype",
                                          'primary_or_metastasis','sex','depmap_public_comments','CCLE_Name')]
    depmap_cellinfo <- depmap_cellinfo[str_detect(depmap_cellinfo$Subtype,'Melanoma'),]
    depmap_cellinfo$primary_or_metastasis[depmap_cellinfo$primary_or_metastasis==""] <- "Uncategorized "
    depmap_cellinfo$sex[depmap_cellinfo$sex==""] <- "Uncategorized "
    depmap_cellinfo$primary_or_metastasis <- factor(depmap_cellinfo$primary_or_metastasis,c("Primary",'Metastasis',"Uncategorized "))
    
    depmap_cellinfo$cell_line_name <- str_replace(depmap_cellinfo$cell_line_name,'-','')
    #depmap_cellinfo$cell_line_name <- str_replace(depmap_cellinfo$cell_line_name,'\\ ','')
    #depmap_cellinfo$cell_line_name <- str_replace(depmap_cellinfo$cell_line_name,'\\.','')
    depmap_cellinfo$cell_line_name <- toupper(depmap_cellinfo$cell_line_name)
    
    drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'-','')
    #drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'\\ ','')
    #drugs_ic50$cell_line_name <- str_replace(drugs_ic50$cell_line_name,'\\.','')
    drugs_ic50$cell_line_name <- toupper(drugs_ic50$cell_line_name)
    
    drugs_ic50 <- drugs_ic50[drugs_ic50$cell_line_name %in% depmap_cellinfo$cell_line_name,] 
    
    #drugs_ic50 <- left_join(drugs_ic50,depmap_cellinfo)
    
    n_distinct(drugs_ic50$Drugs)
    drugs_ic50 <- distinct(drugs_ic50)
    
    drugs_ic50 %>% group_by(Tested_InVitro,Drugs) %>% 
      mutate(Median_IC50 = median(IC50,na.rm = T),
             Std_IC50 = sd(IC50,na.rm = T),
             N_cells = n_distinct(cell_line_name)) ->drugs_ic50
    drugs_ic50$Tested_InVitro <- factor(drugs_ic50$Tested_InVitro,
                                        c("Predicted drug","Approved anti-melanoma","NO-based drug"))
    drugs_ic50$Drugs <- firstup(drugs_ic50$Drugs)
    p12 <- ggplot(drugs_ic50, aes(x=IC50, y=reorder(Drugs,desc(Median_IC50)), 
                                  fill=Tested_InVitro,
                                  color=Tested_InVitro)) + 
      geom_violin(alpha=0.3)+
      geom_point(alpha=0.5)+
      geom_point(aes(x=Median_IC50),color='black')+
      scale_x_continuous(trans='log10',n.breaks = 10)+
      scale_fill_manual(values=c("forestgreen", "#D95F02" ,"#7570B3"))+
      scale_color_manual(values=c("forestgreen", "#D95F02" ,"#7570B3"))+
      ylab("Drugs") +
      #facet_grid(.~primary_or_metastasis)+
      xlab("IC50 in µM") +theme_classic() +
      theme(axis.text.x = element_text(angle = 0,size = 11),
            axis.text.y = element_text(angle = 0,size = 12))+
      geom_vline(xintercept=0.4, linetype="dashed", color = "green", size=0.5)+
      geom_vline(xintercept=1, linetype="dashed", color = "blue", size=0.5)+
      geom_vline(xintercept=10, linetype="dashed", color = "red", size=0.5)+
      ggtitle("IC50 measures in the melanoma cell\nlines across four cell viability databases") +
      guides(fill=guide_legend(title='Drug Type'),color=guide_legend(title='Drug Type'))+
      theme(legend.position = c(0.85,0.62),legend.background = element_blank(),legend.key = element_blank())
    
    p12
    ggsave('Figure/SupplementaryFigureS12_IC50.pdf', p12,units = 'in',width = 10,height = 9,dpi = 300)
    ggsave('Figure/SupplementaryFigureS12_IC50.png', p12,units = 'in',width = 10,height = 9,dpi = 300)
    
    drugs_ic50 %>% group_by(Tested_InVitro,Drugs) %>% 
      summarise(Median_IC50 = median(IC50,na.rm = T),
                Std_IC50 = sd(IC50,na.rm = T),N_cells = n_distinct(cell_line_name)) ->drugs_ic50_summ
    
    write_csv(drugs_ic50,'./Figure/IC50_table.csv')
    write_csv(drugs_ic50_summ,'./Figure/IC50_table_summary.csv')
    library(openxlsx)
    
    write.xlsx(drugs_ic50,'./Figure/IC50_table.xlsx')
    write.xlsx(drugs_ic50_summ,'./Figure/IC50_table_summary.xlsx')