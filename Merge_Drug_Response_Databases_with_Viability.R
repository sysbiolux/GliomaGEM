### Integrating PRISM primary and, Nam et al

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

set_breaks = function(limits) {
  seq(limits[1], limits[2], by = 1)
}

### Analysing SKO drugs with PRISM database for cell lines
# Define glioma cell lines and make a separte column for it
depmap_cellinfo = read.csv("CRISPR/DepMap_22Q1/sample_info_braincancer.csv")
depmap_cellinfo_nonbrain <- depmap_cellinfo[depmap_cellinfo$primary_disease!='Brain Cancer',]
#depmap_cellinfo_nonbrain <- depmap_cellinfo_nonbrain[!depmap_cellinfo_nonbrain$lineage %in% c("peripheral_nervous_system","central_nervous_system"),]
depmap_cellinfo <- depmap_cellinfo[str_detect(depmap_cellinfo$primary_disease,'Brain Cancer'),]

depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID',"lineage_subtype",'lineage_sub_subtype','cell_line_name','primary_or_metastasis','RRID')]

#depmap_cellinfo$Glioma <- ''
#depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma','Glioma'] <- depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma','lineage_sub_subtype']

#depmap_cellinfo <- depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma',] 
depmap_cellinfo <- depmap_cellinfo[depmap_cellinfo$lineage_sub_subtype!='',] 
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-001198','cell_line_name'] <- 'SNB19'
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-001000','cell_line_name'] <- '1321N1'
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-000591','cell_line_name'] <- 'LN235'
depmap_cellinfo$Subtype <- ''
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="astrocytoma" ] <-'AST'
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="glioblastoma" ] <-'GBM'
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="oligodendroglioma" ] <-'ODG'

prism_pr <- read_csv('./Drug_DBs_Resources/PRISM_Repurposing_19Q4/primary-screen-replicate-collapsed-logfold-change.csv')
colnames(prism_pr)[1] <- 'DepMap_ID'
prism_pr <- prism_pr[prism_pr$DepMap_ID %in% depmap_cellinfo$DepMap_ID,]
prism_pr <- as.data.frame(prism_pr)
row.names(prism_pr) <- prism_pr$DepMap_ID

prism_pr_drugs <- read_csv('./Drug_DBs_Resources/PRISM_Repurposing_19Q4/primary-screen-replicate-treatment-info.csv')
prism_pr_drugs <- distinct(prism_pr_drugs[,c('broad_id','name')])

prism_pr %>% pivot_longer(cols = !DepMap_ID,
                            names_to = "broad_conc",values_to = "Log2FoldChange") -> prism_pr_longer

prism_pr_longer %>% separate(broad_conc,sep = '::',into = c('broad_id','conc')) -> prism_pr_longer
prism_pr_longer <- left_join(prism_pr_longer,prism_pr_drugs)

N_celllines = n_distinct(prism_pr_longer$DepMap_ID)
N_drugs <- n_distinct(prism_pr_longer$name)
n_distinct(prism_pr_longer$broad_id)
n_distinct(prism_pr_longer$name)

prism_pr_longer <- prism_pr_longer[,colnames(prism_pr_longer)!='broad_id']
prism_pr_longer <- prism_pr_longer[!is.na(prism_pr_longer$Log2FoldChange),]

prism_pr_longer <- left_join(prism_pr_longer,depmap_cellinfo)

prism_pr_longer$name <- str_replace(prism_pr_longer$name,'\\-',' ')

prism_pr_longer$conc <- lapply(as.numeric(prism_pr_longer$conc) , function(x) signif.ceiling(x, 3))
prism_pr_longer$name_conc <- str_c(prism_pr_longer$name, " (", prism_pr_longer$conc,")")
n_distinct(prism_pr_longer$name_conc)

prism_pr_longer %>% group_by(name_conc)%>% # Subtype
  mutate(#median_logfc = median(Log2FoldChange,na.rm = T),
         FoldChange= (1-(2^Log2FoldChange))*100,
         #median_FoldChange = median(FoldChange,na.rm = T)
         ) -> prism_pr_longer
colnames(prism_pr_longer)
prism_pr_summ <- prism_pr_longer[,c('name','cell_line_name','FoldChange',"Subtype",'conc')]
prism_pr_summ$Database <- 'PRISM Primary'
colnames(prism_pr_summ) <- c('Drugs','Cell_lines',"Reduction_in_Viability",
                             "Subtype",'Dosage_in_µM','Database')
prism_pr_summ$Dosage_in_µM <- as.numeric(prism_pr_summ$Dosage_in_µM)

table(prism_pr_summ$Cell_lines)


##  Nam et al iScience 2021 glioma  drug screening ##############
nam2021 <- readxl::read_xlsx('./Drug_DBs_Resources/Nam_etal_2021_iScience_S_Table1.xlsx')
colnames(nam2021)
colnames(nam2021)[4:6] <- c("Drugs",	"559T",	"592T")
# selecting only glioma cell lines
nam2021 <- nam2021[,4:6]
# lower case
n_distinct(nam2021$Drugs)


nam2021$Drugs <- str_to_lower(nam2021$Drugs)
nam2021 %>% pivot_longer(cols = colnames(nam2021)[2:3],
                         names_to = "Cell_line",values_to = "Reduction_in_Viability",
                         names_ptypes = character()) -> nam2021
nam2021$Subtype <- "GBM"
nam2021$Dosage_in_µM <- 10
nam2021$Database <- "Nam et al 2021"
colnames(nam2021) <- c('Drugs','Cell_lines',"Reduction_in_Viability",
                       "Subtype",'Dosage_in_µM','Database')
nam2021$Drugs <- str_replace(nam2021$Drugs,"vincristine sulfate","vincristine")
nam2021$Drugs <- str_replace(nam2021$Drugs,"2 ,3  - dideoxycytidine","zalcitabine")
nam2021$Drugs <- str_replace(nam2021$Drugs,"5-fluorouracil sulfate","fluorouracil")
nam2021$Drugs <- str_replace(nam2021$Drugs,"vincristine sulfate","vincristine")

Merged_viability <- rbind(prism_pr_summ,nam2021)

glioma_drugs <- read_csv("./Primary_target_databases/Approved_anti_brain_cancers.csv")
queries <- Merged_viability %>% 
  filter(str_detect(Drugs, paste(glioma_drugs$Drugs, collapse = "|")))
unique(queries$Drugs)
Merged_viability$Drugs <- str_replace(Merged_viability$Drugs,"procarbazine hcl","procarbazine")
Merged_viability$Drugs <- str_replace(Merged_viability$Drugs,"doxorubicin hydrochloride","doxorubicin")


write_csv(Merged_viability,'./Merged_Viability_databases.csv')
# 
# n_distinct(nam2021$Cell_line)
# 
# # take average [ median not going to work with two]
# # Calculate the n percentile of the median rank
# nam2021 %>% group_by(Drugs) %>% 
#   #filter(is.na(TYPE2)) %>%
#   mutate(median_per_drug=median(Viability))->nam2021
# 
# nam2021_sko <- left_join(nam2021,sko_df_)
# 
# nam2021_sko <- distinct(nam2021_sko[,c('Drugs','median_per_drug','TYPE2')])
# nam2021_sko$Nam2021_etal<- ntile(-nam2021_sko$median_per_drug,100)
# 
# Ranking_df_5 <- distinct(nam2021_sko[!is.na(nam2021_sko$TYPE2),])
# nam2021_sko
# ###### Merge all HTS 
# library(plyr)
# colnames(Ranking_df_1)[1] <- 'Drugs' 
# colnames(Ranking_df_2)[1] <- 'Drugs' 
# 
# #ranking_all <- join_all(list(sko_df,Ranking_df_1,Ranking_df_2,Ranking_df_3,Ranking_df_4[,c('Drugs','Bell_etal_Folllowup')],
# #                         Ranking_df_5[,c('Drugs','Nam2021_etal')]),by='Drugs',type='left')
# ranking_all <- join_all(list(sko_df_[,c('Drugs','TYPE2')],Ranking_df_1,Ranking_df_2, Ranking_df_5[,c('Drugs','Nam2021_etal')]),by='Drugs',type='left')
# ranking_all %>% distinct() %>%
#   pivot_longer(cols = colnames(ranking_all)[3:5],
#                names_to = "HTS_database",values_to = "Rank",
#                names_ptypes = character()) -> ranking_all
# ranking_all <- ranking_all[!is.na(ranking_all$Rank),]
# ranking_all %>%  group_by(Drugs,TYPE2,HTS_database) %>% slice(which.min(Rank)) ->ranking_all
# 
# 
# detach("package:plyr", unload = TRUE)
# 
# ranking_all %>% group_by(Drugs) %>% mutate(median_rank = median(Rank)) -> ranking_all
# ranking_all$Predicted <- 'Yes'
# ranking_all$Predicted[ranking_all$TYPE2=="Approved Anti-brain cancer"] <- 'Approved Anti\nbrain cancer'
# 
# ranking_all$TYPE2[ranking_all$TYPE2=="Approved Anti-brain cancer"] <- 'Single Drug Deletion'
# ranking_all$TYPE2 <- factor(ranking_all$TYPE2,c("Single Drug Deletion","Double Drug Deletion"))
# 
# ggplot(ranking_all,aes(y=fct_reorder(Drugs,desc(median_rank)),alpha=Predicted,
#                        x = Rank,
#                        fill=HTS_database))+
#   geom_bar(stat = "identity", position=position_dodge())+
#   #geom_text(aes(label=Rank),size=3.5)+
#   ggtitle("Rank of the predicted drugs in 3 in vitro high-throughput drug screening databases") +
#   xlab("Rank in nth percentile among all screened drugs") +
#   ylab("Drugs") +theme_classic() +
#   scale_alpha_discrete(range=c(0.5,1),limits=rev)+
#   facet_wrap(.~TYPE2,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12))

# prism_pr_sko <- left_join(sko_df_[,c('Drugs','Subtype','TYPE2')], prism_pr_drugs,by=c('Drugs'='name'))
# unique(prism_pr_sko$Drugs[is.na(prism_pr_sko$broad_id)]) # drugs not in the database
# prism_pr_drugs_conc <- as.data.frame(colnames(prism_pr))
# prism_pr_drugs_conc %>% separate(`colnames(prism_pr)`,sep = '::',into = c('broad_id','conc')) ->prism_pr_drugs_conc
# prism_pr_sko <- left_join(prism_pr_sko,prism_pr_drugs_conc)
# prism_pr_sko$broad_conc <- str_c(prism_pr_sko$broad_id,'::',prism_pr_sko$conc,'::HTS') 
# prism_pr_sko_ <- prism_pr[,colnames(prism_pr)  %in% prism_pr_sko$broad_conc]
# #prism_pr_sko <- prism_pr_sko[!is.na(prism_pr_sko$broad_id),]
# 
# prism_pr_sko_$DepMap_ID  <- rownames(prism_pr_sko_)
# prism_pr_sko_ %>% pivot_longer(cols = !DepMap_ID,
#                                names_to = "broad_conc",values_to = "Log2FoldChange") -> prism_pr_sko_
# 
# prism_pr_sko <- left_join(depmap_cellinfo,prism_pr_sko)
# prism_pr_sko_final <- left_join(prism_pr_sko,prism_pr_sko_) #Joining, by = c("DepMap_ID", "broad_conc")
# 
# ## The number of cell line used for each drug
# prism_pr_sko_final %>% group_by(Subtype,Drugs) %>% 
#   summarise(n_celllines = n_distinct(DepMap_ID[!is.na(Log2FoldChange)])) ->prism_pr_sko_summ
# 
# # Remove ODG since it has zero cell lines for all drugs
# prism_pr_sko_summ <- prism_pr_sko_summ[prism_pr_sko_summ$Subtype!='ODG',]
# 
# ggplot(prism_pr_sko_summ,aes(y=reorder_within(Drugs,n_celllines,Subtype),x = as.integer(n_celllines),fill=Subtype))+
#   geom_bar(stat = "identity", position=position_dodge())+
#   ggtitle("Number of cell lines tested for SKO drugs in Primary PRISM database") +
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
#   scale_y_reordered()+
#   xlab("Number of tested cell lines") +
#   ylab("Drugs") +theme_classic() +
#   facet_wrap(.~Subtype,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12))
# 
# 
# signif.ceiling(1.00000235,3)
# # Ceiling the conc to the 2nd decimal
# prism_pr_sko_final$conc <- lapply(as.numeric(prism_pr_sko_final$conc) , function(x) signif.ceiling(x, 3))
# 
# # Remove drugs with no log FC
# prism_pr_sko_final <- prism_pr_sko_final[!is.na(prism_pr_sko_final$Log2FoldChange),]
# 
# # Replace broad_id with names + conc
# prism_pr_sko_final$Drugs_conc <- str_c(prism_pr_sko_final$Drugs,'_',prism_pr_sko_final$conc)
# 
# 
# prism_pr_sko_final %>%
#   dplyr::group_by(Drugs_conc, cell_line_name,Subtype,TYPE2) %>%
#   dplyr::summarise(Log2FoldChange = mean(Log2FoldChange)) %>%
#   group_by(Drugs_conc,Subtype) %>% mutate(average_per_drug = mean(Log2FoldChange)) %>%
#   group_by(cell_line_name,Subtype) %>% mutate(average_per_cellline = mean(Log2FoldChange)) ->prism_pr_sko_final_
# 
# # Ceiling the conc to the 2nd decimal
# prism_pr_sko_final_$Log2FoldChange_ <- lapply(as.numeric(prism_pr_sko_final_$Log2FoldChange) , function(x) signif.ceiling(x, 3))
# 
# min_value <- min(prism_pr_sko_final_$Log2FoldChange)
# max_value <- max(prism_pr_sko_final_$Log2FoldChange)
# p <- ggplot(prism_pr_sko_final_,aes(y=reorder_within(Drugs_conc,desc(average_per_drug),Subtype),
#                                     x = reorder_within(cell_line_name,average_per_cellline,Subtype),fill=Log2FoldChange))+
#   geom_tile()+#stat = "identity", position=position_dodge())+
#   ggtitle("Drug sensitivity of the predicted drugs in Primary PRISM database") +
#   geom_text(aes(label=Log2FoldChange_),size=2)+
#   #scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
#   #scale_fill_brewer(palette="Dark2")+
#   #scale_fill_distiller(palette = "OrRd", direction = 1,trans ='reverse')+
#   scale_fill_gradientn(colours = c("cyan", "white", "red"), #na.value = na.value.forplot,
#                        values = scales::rescale(c(max_value, 0, min_value)))+
#   #scale_fill_continuous(na.value="grey")+
#   scale_x_reordered()+
#   scale_y_reordered()+
#   xlab("Tested cell lines") +
#   ylab("Drugs with conc in micro M") +theme_classic() +
#   facet_wrap(TYPE2~Subtype,scales = 'free')+ # ,space="free"
#   #  scale_x_discrete(expand = c(0, 0.5)) + # change additive expansion from default 0.6 to 0.5
#   
#   theme(axis.text.x = element_text(angle = 90,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12),
#         panel.background = element_rect(fill = "grey92", colour = NA))
# p
# 
# ggsave('Figs/PRISM_1st_drugs.png', p,units = 'in',width = 20,height = 14,dpi = 300)
# 
# 
# ## Define the rank of the sko drugs in the glioma cell lines (overall rank not by subtype)
# prism_pr %>% pivot_longer(cols = !DepMap_ID,names_to = "broad_conc",values_to = "Log2FoldChange") -> prism_pr_longer
# prism_pr_longer %>% separate(broad_conc,sep = '::',into = c('broad_id','conc')) -> prism_pr_longer
# prism_pr_longer <- left_join(prism_pr_longer,prism_pr_drugs)
# 
# prism_pr_longer <- left_join(prism_pr_longer,sko_df_,by=c('name'='Drugs')) 
# 
# # selecting only glioma cell lines
# prism_pr_longer <- prism_pr_longer[prism_pr_longer$DepMap_ID%in% depmap_cellinfo$DepMap_ID,]
# #prism_pr_longer <- left_join(prism_pr_longer,depmap_cellinfo[,c('DepMap_ID','Subtype')])
# 
# 
# length(unique(prism_pr_longer$DepMap_ID))
# 
# prism_pr_longer %>% group_by(name)%>% # Subtype
#   mutate(median_logfc = median(Log2FoldChange,na.rm = T)) %>%
#   group_by(DepMap_ID)%>%
#   dplyr::arrange(median_logfc) %>%
#   rank_in_group(.)->prism_pr_longer_rank
# 
# # Calculate the n percentile of the median rank
# prism_pr_longer_rank$PRISM_Primary <- ntile(prism_pr_longer_rank$median_logfc,100)
# Ranking_df_1 <- distinct(prism_pr_longer_rank[!is.na(prism_pr_longer_rank$TYPE2),c('name','PRISM_Primary')])
# 
# length(unique(prism_pr_longer_rank$name))
# prism_pr_longer_rank %>% group_by(name) %>%# 
#   mutate(rank=max(rank))%>% distinct() ->prism_pr_longer_rank
# 
# prism_pr_longer_rank <- distinct(prism_pr_longer_rank[,c('name','Subtype','TYPE2','median_logfc','rank')])
# 
# prism_pr_longer_rank <- na.omit(prism_pr_longer_rank)
# 
# p_ <- ggplot(prism_pr_longer_rank,aes(y=reorder(name,desc(rank)),
#                                       #fill=Subtype,#
#                                       x = median_logfc))+
#   geom_bar(stat = "identity", position=position_dodge())+
#   geom_text(aes(label=rank),size=4.5)+
#   ggtitle("Rank of the predictd drugs in the glioma cell lines\nout of 4518 drugs in the Primary PRISM database") +
#   #scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
#   #scale_fill_brewer(palette="Dark2")+
#   #scale_y_reordered()+
#   xlab("Median Log2 Fold Change") +
#   ylab("Drugs") +theme_classic() +
#   #facet_wrap(.~,scales = 'free')+
#   facet_wrap(.~TYPE2,scales = 'free')+
#   
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 12))
# p_
# 
# ggsave('Figs/PRISM_1st_drugs_rank.png', p_,units = 'in',width = 20,height = 14,dpi = 300)

