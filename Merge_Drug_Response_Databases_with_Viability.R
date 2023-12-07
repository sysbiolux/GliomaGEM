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
