### Integrating 3 xenograft databases for GBM
## Bell et al initial
## Bell et al followup
## Stathias et al 2018, dowloaded from https://data.mendeley.com/datasets/yz8m28gj6r/1

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


### Stathias et al 2018
stathias2018 <- readxl::read_xlsx('Drug_DBs_Resources/Xenograft_data/Stathias_etal_2018/Supplementary File 1.xlsx',
                                  sheet = "FDA approved screen")
colnames(stathias2018)[2:8] <- stathias2018[1,2:8]
stathias2018 <- stathias2018[2:nrow(stathias2018),1:8]
stathias2018 %>% pivot_longer(!Compound,names_to = "PDX",values_to = "Viability_Reduction") ->stathias2018
stathias2018$Database <- "Stathias et al 2018"
colnames(stathias2018)[1] <- "Drugs"
stathias2018$Dosage_in_µM <- 1

### Find the drug sensitivity from Bell et al 2018 https://doi.org/10.1158/1541-7786.MCR-17-0397
bell2018 <- readxl::read_xlsx('Drug_DBs_Resources//Xenograft_data/Bell_etal_Mol_Cancer_Res_2018_Supple_Table1.xlsx')
bell2018_2 <- readxl::read_xlsx('Drug_DBs_Resources/Xenograft_data/Bell_etal_Mol_Cancer_Res_2018_Supple_Table2.xlsx')
bell2018_2_dico <- bell2018_2[4:nrow(bell2018_2),c(1,5)]
colnames(bell2018_2_dico)[1] <- 'SBI_ID'
colnames(bell2018_2_dico)[2] <- 'Drugs_2'

colnames(bell2018)[1] <- 'SBI_ID'
colnames(bell2018)[2] <- 'Drugs'
colnames(bell2018)[3:7] <- bell2018[2,3:7]
colnames(bell2018)[8:9] <- bell2018[1,8:9]
bell2018 <- bell2018[3:nrow(bell2018),]

bell2018 %>% pivot_longer(cols = c("GBM102","GBM84", "GBM64","GBM46" , "GBM12"  , "U87", "HFF") ,
                          names_to = "Cell_line",values_to = "Viability") -> bell2018
bell2018 <- bell2018[bell2018$Cell_line!='HFF',]
intersect(bell2018$SBI_ID,bell2018_2_dico$SBI_ID)
bell2018 <- left_join(bell2018,bell2018_2_dico)
bell2018$Drugs[is.na(bell2018$Drugs)] <- bell2018$Drugs_2[is.na(bell2018$Drugs)]

#bell2018$Drugs[is.na(bell2018$Drugs)] <- bell2018$`Compound Name`[is.na(bell2018$Drugs)] 

bell2018$Drugs  <- as.character(bell2018$Drugs )
bell2018$Viability  <- as.numeric(bell2018$Viability )

bell2018$Drugs <- str_replace(str_to_lower(bell2018$Drugs),' ','-')
bell2018$Drugs <- str_replace(bell2018$Drugs,"- (5-fu)","")
bell2018$Drugs <- str_replace(bell2018$Drugs,"-hydrochloride","")
#gbm_drugs <- sko_df[sko_df$Subtype=='GBM','Drugs']
bell2018_summ <- bell2018[,c('Drugs','Cell_line','Viability')]
bell2018_summ$Viability <- 100- bell2018_summ$Viability*100
colnames(bell2018_summ) <- c('Drugs','PDX','Viability_Reduction')
bell2018_summ$Database <- "Bell et al 2018, initial"
bell2018_summ$Dosage_in_µM <- 10
## Followup screening
### Find the drug sensitivity from Bell et al 2018 https://doi.org/10.1158/1541-7786.MCR-17-0397

#bell2018_2 <- bell2018_2_df
colnames(bell2018_2)[6:53] <- bell2018_2[2,6:53]
colnames(bell2018_2)[seq(6,52,9)+1] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+2] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+3] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+4] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+5] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+6] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+7] <- colnames(bell2018_2)[seq(6,52,9)]
colnames(bell2018_2)[seq(6,52,9)+8] <- colnames(bell2018_2)[seq(6,52,9)]


colnames(bell2018_2)[seq(6,50,3)] <- str_c(colnames(bell2018_2)[seq(6,50,3)] ,'_',bell2018_2[3, seq(6,50,3)])
colnames(bell2018_2)[seq(6,50,3)+1] <- str_c(colnames(bell2018_2)[seq(6,50,3)],'_2')
colnames(bell2018_2)[seq(6,50,3)+2] <- str_c(colnames(bell2018_2)[seq(6,50,3)],'_3')


bell2018_2 <- bell2018_2[4:nrow(bell2018_2), 5:50]

colnames(bell2018_2)[1] <- 'Drugs'
bell2018_2[,2:46] <- lapply(bell2018_2[,2:46],as.numeric)
bell2018_2 %>% pivot_longer(cols = colnames(bell2018_2)[2:46],
                            names_to = "Cell_line",values_to = "Viability",
                            names_ptypes = character()) -> bell2018_2

bell2018_2$Drugs <- str_to_lower(bell2018_2$Drugs)
bell2018_2$Drugs <- str_replace_all(bell2018_2$Drugs," ","-")

bell2018_2$Drugs[bell2018_2$Drugs=="fluorouracil--(5-fu)"] <- "fluorouracil"
bell2018_2$Drugs[bell2018_2$Drugs=="lomustine-(ccnu)"] <- "lomustine"

bell2018_2$Drugs <- str_replace(bell2018_2$Drugs,"-hcl","")
bell2018_2$Drugs <- str_replace(bell2018_2$Drugs,"trametinib-(gsk1120212)","trametinib")

bell2018_2_summ <- bell2018_2[,c('Drugs','Cell_line','Viability')]
bell2018_2_summ$Viability <- 100- bell2018_2_summ$Viability 
colnames(bell2018_2_summ) <- c('Drugs','PDX','Viability_Reduction')
bell2018_2_summ$Database <- "Bell et al 2018, followup"
bell2018_2_summ$Dosage_in_µM <- str_split(bell2018_2_summ$PDX,'_',simplify = TRUE)[,2]
bell2018_2_summ$Dosage_in_µM <- str_split(bell2018_2_summ$Dosage_in_µM,' ',simplify = TRUE)[,1]
bell2018_2_summ$PDX <- str_split(bell2018_2_summ$PDX,'_',simplify = TRUE)[,1]
bell2018_2_summ %>% group_by(Drugs,PDX, Dosage_in_µM,Database) %>% 
  summarise(Viability_Reduction = median(Viability_Reduction)) -> bell2018_2_summ

### Merge the three xenograft databases 
colnames <- colnames(stathias2018)
integrated_xeno <- rbind(stathias2018,bell2018_summ[,colnames],bell2018_2_summ[,colnames])
integrated_xeno$Viability_Reduction <- as.numeric(integrated_xeno$Viability_Reduction )
integrated_xeno$Dosage_in_µM <- as.numeric(integrated_xeno$Dosage_in_µM )

integrated_xeno$Drugs <- str_to_lower(integrated_xeno$Drugs )
integrated_xeno$Drugs[integrated_xeno$Drugs=="cyclophosphamide monohydrate"] <- "cyclophosphamide"
integrated_xeno$Drugs[integrated_xeno$Drugs=="everolimus (rad001)"] <- "everolimus"
integrated_xeno$Drugs[integrated_xeno$Drugs=="doxorubicin (adriamycin)"] <- "doxorubicin"
integrated_xeno$Drugs[integrated_xeno$Drugs=="fludarabine (fludara)"] <- "fludarabine-phosphate"
integrated_xeno$Drugs[integrated_xeno$Drugs=="fludarabine phosphate (fludara)"] <- "fludarabine-phosphate"
integrated_xeno$Drugs[integrated_xeno$Drugs=="adrucil (fluorouracil)"] <- "fluorouracil"
integrated_xeno$Drugs[integrated_xeno$Drugs=="gemcitabine hcl (gemzar)"] <- "gemcitabine"
integrated_xeno$Drugs[integrated_xeno$Drugs=="lomustine (ceenu)"] <- "lomustine"
integrated_xeno$Drugs[integrated_xeno$Drugs=="dorzolamide hcl"] <- "dorzolamide"
integrated_xeno$Drugs[integrated_xeno$Drugs=="vincristine-sulfate"] <- "vincristine"
integrated_xeno$Drugs[integrated_xeno$Drugs=="trametinib-(gsk1120212)"] <- "trametinib"

integrated_xeno$Drugs<- str_split(integrated_xeno$Drugs," ",simplify = T)[,1]
write_csv(integrated_xeno,"./Drug_DBs_Resources/Integrated_xenografts_data.csv")

# ##############
# bell2018 %>% group_by(Drugs) %>% 
#   mutate(median_viability  = median(Viability[Cell_line!='U87'])) %>%
#   group_by(Cell_line)%>%
#   dplyr::arrange(median_viability) %>%
#   rank_in_group(.)->bell2018
# 
# unique(bell2018$Drugs[str_detect(bell2018$Drugs , paste0(sko_df_$Drugs,collapse = '|'))])
# bell2018_ <- distinct(bell2018[,c('Drugs','median_viability')])
# bell2018_$Bell_etal_Initial<- ntile(bell2018_$median_viability,100)
# bell2018_ <- left_join(bell2018_,sko_df_)
# Ranking_df_3 <- distinct(bell2018_[!is.na(bell2018_$TYPE2),c('Drugs','Bell_etal_Initial')])
# 
# #bell2018_sko <- bell2018[bell2018$Drugs %in% gbm_drugs$Drugs,]
# bell2018_sko <- left_join(bell2018,sko_df_)
# #x <- bell2018_sko[!is.na(bell2018_sko$TYPE2),]
# #bell2018_sko <- bell2018[bell2018$Drugs %in% union(gbm_drugs$Drugs,"gemcitabine-hydrochloride"),]
# 
# # Ceiling the measure to the 2nd decimal
# bell2018_sko$Measure_value_ <- lapply(as.numeric(bell2018_sko$Viability) , function(x) signif.ceiling(x, 3))
# 
# bell2018_sko$Celltype <- 3#'Cell line'
# bell2018_sko$Celltype[bell2018_sko$Cell_line %in%c("GBM102","GBM84", "GBM64") ] <- 1#'Mesenchymal PDX'
# bell2018_sko$Celltype[bell2018_sko$Cell_line %in%c("GBM46" , "GBM12") ] <- 2#'Proneural PDX'
# 
# bell2018_sko <- bell2018_sko[!is.na(bell2018_sko$TYPE2),]
# 
# 
# ggplot(bell2018_sko[bell2018_sko$TYPE2 %in% c('Single Drug Deletion','Approved Anti-brain cancer'),],aes(y=fct_reorder(Drugs,desc(median_viability)),
#                                                                                                          x = fct_reorder(Cell_line,Celltype),fill=Viability))+
#   geom_tile()+
#   ggtitle("Drug sensitivity of 650 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   geom_text(aes(label=Measure_value_),size=4)+
#   scale_fill_distiller(palette = "OrRd", direction = -1)+
#   labs(fill="Viability to control (%)")+
#   xlab("Tested GBM PDX models / cell lines") +
#   ylab("Drugs") +theme_classic() +
#   facet_wrap(TYPE2~.,scales = 'free',nrow = 2)+
#   theme(axis.text.x = element_text(angle = 0,size = 13),
#         axis.text.y = element_text(angle = 0,size = 13),
#         panel.background = element_rect(fill = "grey92", colour = NA))
# 
# ggplot(bell2018_sko[bell2018_sko$TYPE2=='Double Drug Deletion',],aes(y=fct_reorder(Drugs,desc(median_viability)),
#                                                                      x = fct_reorder(Cell_line,Celltype),fill=Viability))+
#   geom_tile()+
#   ggtitle("Drug sensitivity of 650 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   geom_text(aes(label=Measure_value_),size=4)+
#   scale_fill_distiller(palette = "OrRd", direction = -1)+
#   labs(fill="Viability to control (%)")+
#   xlab("Tested GBM PDX models / cell lines") +
#   ylab("Drugs") +theme_classic() +
#   facet_wrap(TYPE2~.,scales = 'free',nrow = 2)+
#   theme(axis.text.x = element_text(angle = 0,size = 13),
#         axis.text.y = element_text(angle = 0,size = 13),
#         panel.background = element_rect(fill = "grey92", colour = NA))
# 
# ### Follow up Bell et al
# unique(bell2018_2$Drugs[str_detect(bell2018_2$Drugs , paste0(sko_df_$Drugs,collapse = '|'))])
# 
# 
# bell2018_2_sko <- left_join(bell2018_2,sko_df_)
# bell2018_2_sko$Cellline <- str_split(bell2018_2_sko$Cell_line,'_',simplify = T)[,1]
# bell2018_2_sko$Conc <- str_split(bell2018_2_sko$Cell_line,'_',simplify = T)[,2]
# bell2018_2_sko$Replicate <- str_split(bell2018_2_sko$Cell_line,'_',simplify = T)[,3]
# 
# 
# # Calculate the n percentile of the median rank
# bell2018_2_sko %>% group_by(Drugs) %>% 
#   #filter(is.na(TYPE2)) %>%
#   mutate(median_per_drug=median(Viability))->bell2018_2_sko_
# bell2018_2_sko_ <- distinct(bell2018_2_sko_[,c('Drugs','median_per_drug','TYPE2')])
# bell2018_2_sko_$Bell_etal_Folllowup<- ntile(bell2018_2_sko_$median_per_drug,100)
# Ranking_df_4 <- distinct(bell2018_2_sko_[!is.na(bell2018_2_sko_$TYPE2),])
# 
# 
# bell2018_2_sko %>% group_by(Drugs,Conc) %>% mutate(median_viability=median(Viability))  %>% 
#   group_by(Drugs,Conc,Cellline) %>% mutate(mean_per_replicate=mean(Viability))->bell2018_2_sko
# 
# 
# # Ceiling the measure to the 2nd decimal
# bell2018_2_sko$Measure_value_ <- lapply(as.numeric(bell2018_2_sko$Viability) , function(x) signif.ceiling(x, 3))
# bell2018_2_sko$mean_per_replicate_ <- lapply(as.numeric(bell2018_2_sko$mean_per_replicate) , function(x) signif.ceiling(x, 3))
# 
# bell2018_2_sko$Celltype <-NA
# bell2018_2_sko$Celltype[bell2018_2_sko$Cellline %in%c("GBM91","GBM102") ] <- 'Mesenchymal'#'Mesenchymal PDX'
# bell2018_2_sko$Celltype[bell2018_2_sko$Cellline %in%c("SF7300","GBM59"	,	"GBM116") ] <- 'Proneural'#'Proneural PDX'
# bell2018_2_sko$Cellline <- factor(bell2018_2_sko$Cellline,c("GBM91","GBM102","SF7300","GBM59","GBM116"))
# 
# bell2018_2_sko$Drugs_Conc_Rep <- str_c(bell2018_2_sko$Drugs,'_',bell2018_2_sko$Conc,'_',bell2018_2_sko$Replicate)
# bell2018_2_sko$Drugs_Conc <- str_c(bell2018_2_sko$Drugs,'_',bell2018_2_sko$Conc)
# 
# bell2018_2_sko <- bell2018_2_sko[!is.na(bell2018_2_sko$TYPE2),]
# bell2018_2_sko$Drugs_Conc <- str_replace(bell2018_2_sko$Drugs_Conc,'mM','µM')
# # ggplot(bell2018_2_sko,aes(y=fct_reorder(Drugs_Conc,desc(median_viability)),
# #                         x = fct_reorder(Cellline,Celltype),
# #                         fill=mean_per_replicate))+
# #   geom_tile()+
# #   ggtitle("Drug sensitivity of 120 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
# #   #guides(fill=guide_legend(title=str_to_upper(MEASURE),nrow = 5))+
# #   geom_text(aes(label=mean_per_replicate_),size=3)+
# #   scale_fill_distiller(palette = "OrRd", direction = -1)+
# #   labs(fill='Average viability % \nto control')+
# #   #  scale_x_reordered()+
# #   #  scale_y_reordered()+
# #   xlab("Tested GBM PDX models") +
# #   ylab("Drugs") +theme_classic() +
# #   facet_wrap(TYPE2~.,scales = 'free',nrow = 2)+
# #   theme(axis.text.x = element_text(angle = 0,size = 11),
# #         axis.text.y = element_text(angle = 0,size = 12),
# #         panel.background = element_rect(fill = "grey92", colour = NA))
# 
# ### Dose analysis of Bell followup
# library(ggpmisc)
# 
# bell2018_2_sko$Conc <- as.numeric(str_replace(bell2018_2_sko$Conc, ' mM',''))
# bell2018_2_sko <- distinct(bell2018_2_sko[,colnames(bell2018_2_sko)!='Subtype'])
# bell2018_2_sko <- distinct(bell2018_2_sko[,colnames(bell2018_2_sko)!='DelRxns'])
# bell2018_2_sko <- distinct(bell2018_2_sko[,colnames(bell2018_2_sko)!='SYMBOL'])
# 
# 
# bell2018_2_sko$Group <- str_c(bell2018_2_sko$Drugs,' (',bell2018_2_sko$Cellline,')')
# 
# b_split = split(bell2018_2_sko$Viability,bell2018_2_sko$Group)
# c_split = split(log10(bell2018_2_sko$Conc),bell2018_2_sko$Group)
# 
# getR2 = function(ip,ip2){
#   model = lm(ip~ip2)
#   return(summary(model)$adj.r.squared)
# }
# 
# r2_df = as.data.frame(mapply(getR2,b_split,c_split))
# r2_df$Group <- rownames(r2_df)
# colnames(r2_df)[1] <- 'R2'
# bell2018_2_sko <- left_join(bell2018_2_sko,r2_df)
# 
# r2_df$Drugs <- str_split(r2_df$Group," \\(",simplify = T)[,1]
# r2_df$Cell_line <- str_split(r2_df$Group," \\(",simplify = T)[,2]
# r2_df$Cell_line <- str_replace(r2_df$Cell_line,"\\)","")
# r2_df %>%group_by(Drugs)%>% mutate(median_R2 = median(R2)) ->r2_df
# unique(bell2018_2_sko$Drugs )
# bell2018_2_sko$Drugs <-factor(bell2018_2_sko$Drugs,unique(r2_df$Drugs[order(r2_df$median_R2,decreasing = T)]))
# bell2018_2_sko$Drugs <-factor(bell2018_2_sko$Drugs,c("gemcitabine","clofarabine","cladribine","mercaptopurine","melphalan","arsenic-trioxide","fludarabine-phosphate","decitabine","fluorouracil","celecoxib"))
# ggplot(bell2018_2_sko[bell2018_2_sko$TYPE2=='Single Drug Deletion',],aes(y=Viability,x = Conc,color=Drugs))+#
#   geom_point()+
#   ggtitle("Dose-viability relationship among the 120 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   stat_poly_line() +
#   stat_poly_eq(label.x = "right",label.y = "top")+
#   stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-val = ", signif(..p.value.., digits = 2), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
#   scale_x_continuous(trans = 'log10')+
#   xlab("Dosage") +
#   ylab("Viability") +theme_classic() +
#   facet_grid(Celltype+Cellline~Drugs,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 10),
#         axis.text.y = element_text(angle = 0,size = 12),
#         panel.background = element_rect(fill = "grey92", colour = NA))+
#   scale_color_discrete(guide="none")
# 
# bell2018_2_sko$Drugs <- str_split(bell2018_2_sko$Drugs_Conc,'_',simplify = T)[,1]
# ggplot(bell2018_2_sko[bell2018_2_sko$TYPE2=='Approved Anti-brain cancer',],aes(y=Viability,x = Conc,color=Drugs))+#
#   geom_point()+
#   ggtitle("Dose-viability relationship among the 120 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   stat_poly_line() +
#   stat_poly_eq(label.x = "right",label.y = "top")+
#   stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-val = ", signif(..p.value.., digits = 2), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
#   scale_x_continuous(trans = 'log10')+
#   xlab("Dosage") +
#   ylab("Viability") +theme_classic() +
#   facet_grid(Celltype+Cellline~Drugs,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 10),
#         axis.text.y = element_text(angle = 0,size = 12),
#         panel.background = element_rect(fill = "grey92", colour = NA))+
#   scale_color_discrete(guide="none")
# 
# ggplot(bell2018_2_sko[bell2018_2_sko$TYPE2=='Double Drug Deletion',],aes(y=Viability,x = Conc,color=Drugs))+#
#   geom_point()+
#   ggtitle("Dose-viability relationship among the 120 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   stat_poly_line() +
#   stat_poly_eq(label.x = "right",label.y = "top")+
#   stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-val = ", signif(..p.value.., digits = 2), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
#   scale_x_continuous(trans = 'log10')+
#   xlab("Dosage") +
#   ylab("Viability") +theme_classic() +
#   facet_grid(Celltype+Cellline~Drugs,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 10),
#         axis.text.y = element_text(angle = 0,size = 12),
#         panel.background = element_rect(fill = "grey92", colour = NA))+
#   scale_color_discrete(guide="none")
# 
# colnames(bell2018_2_sko)
# bell2018_2_sko_summ <- distinct(bell2018_2_sko[,c("Drugs","Cellline","TYPE2","R2","Celltype")])
# #bell2018_2_sko_summ <- bell2018_2_sko_summ[bell2018_2_sko_summ$TYPE2=='Single Drug Deletion',]
# 
# drugs_odgonly <- setdiff(sko_df$Drugs[sko_df$Subtype=='ODG'],sko_df$Drugs[sko_df$Subtype!='ODG'])
# allglioma_drugs <- intersect(sko_df$Drugs[sko_df$Subtype=='AST'],sko_df$Drugs[sko_df$Subtype=='GBM'])
# allglioma_drugs <- intersect(allglioma_drugs,sko_df$Drugs[sko_df$Subtype=='ODG'])
# bell2018_2_sko_summ$Prediction <- ''
# bell2018_2_sko_summ$Prediction[bell2018_2_sko_summ$Drugs %in% drugs_odgonly] <-'Only for ODG'
# bell2018_2_sko_summ$Prediction[bell2018_2_sko_summ$Drugs %in% allglioma_drugs] <-'All Subtypes'
# bell2018_2_sko_summ %>%group_by(Drugs)%>% mutate(median_R2 = median(R2)) ->bell2018_2_sko_summ
# 
# 
# ggplot(bell2018_2_sko_summ,aes(x=R2,y = reorder(Drugs,median_R2),color=Celltype))+#
#   geom_point()+
#   ggtitle("Dose-viability relationship among the 120 drugs in GBM PDX models\nfrom Bell et al 2018 Mol Cancer Res") +
#   #stat_poly_line() +
#   #stat_poly_eq(label.x = "right") +
#   #stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-val = ", signif(..p.value.., digits = 2), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
#   #scale_x_continuous(trans = 'log10')+
#   xlab("R2 of dose-viability") +
#   ylab("Drugs") +theme_classic() +
#   facet_wrap(.~TYPE2,scales = 'free')+
#   theme(axis.text.x = element_text(angle = 0,size = 10),
#         axis.text.y = element_text(angle = 0,size = 12),
#         panel.background = element_rect(fill = "grey92", colour = NA))
# #scale_color_discrete(guide="none")
# 


#bell_df$Xenograft_or_Animal <- 'Xenograft'
#bell_df$Effect_with_TMZ <- 'Not tested'


