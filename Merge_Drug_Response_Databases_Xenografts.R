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
