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

library(qdapRegex)
library(tidyverse)
setwd("./GliomGEM")
## File has replaced " ! " with "_" manually
df <- get_OBO('CRISPR/DepMap_22Q1/cellosaurus.obo.txt', extract_tags = "everything")
df2 <-   simplify2array(df) 
df2 <- as.data.frame(df2)
df2[7467,'xref']
df2[df2$id =='CVCL_C3IQ','xref']

df2 %>%
  as_tibble() %>% #Convert to tibble
  #select(id,comment) %>% #select HPO ID and xref
  unnest(c(id,comment)) %>% #unnest list columns
  separate(comment, into = c("Ontology","Term"), sep = "\\:") %>% #separate ontology from code
  pivot_wider(id_cols = id, names_from = "Ontology",
              values_from = Term,
              values_fn = \(x)paste(x,collapse = ";")) -> df2_comment
df2 %>% #Convert to array
  as_tibble() %>% #Convert to tibble
  select(id,xref) %>% #select HPO ID and xref
  unnest(c(id,xref))  %>% #unnest list columns
  separate(xref, into = c("Ontology","Term"), sep = ":") %>%   
  pivot_wider(id_cols = id, names_from = "Ontology",
              values_from = Term,
              values_fn = \(x)paste(x,collapse = ";")) -> df2_xref

unique(df2_xref$NCIt)[1000:2000]
subtypes <- c('glioblastoma','astrocytoma','oligodendroglioma')
df2_xref$NCIt <- str_to_lower(df2_xref$NCIt)
colnames(df2_xref)

df4 <- df2_xref[str_detect(df2_xref$NCIt,paste0(subtypes,collapse='|')),]
df4 %>% drop_na(id) -> df4
df4$NCIt
unique(df4$NCBI_TaxID)
df4 <- df4[str_detect(df4$NCBI_TaxID,"9606_Homo sapiens"),]

df4$TYPE <- ''
df4$TYPE[str_detect(df4$NCIt,'astrocytoma')] <- 'Astrocytoma'
df4$TYPE[str_detect(df4$NCIt,'glioblastoma')] <- 'Glioblastoma'
df4$TYPE[str_detect(df4$NCIt,'oligodendroglioma')] <- 'Oligodendroglioma'

df2_comment$id <- as.character(df2_comment$id)
df4 <- left_join(df4[,c('id','TYPE')],df2_comment)
df4 <- df4[,c('id','TYPE')]
##############

colnames(df4) <- c('RRID','Subtype2')
df4$Subtype_cellosaurus <- ""
df4$Subtype_cellosaurus[str_detect(df4$Subtype2,'Astrocytoma')] <- 'AST'
df4$Subtype_cellosaurus[str_detect(df4$Subtype2,'Glioblastoma')] <- 'GBM'
df4$Subtype_cellosaurus[str_detect(df4$Subtype2,'Oligodendroglioma')] <- 'ODG'


### Using DepMpa metadata
# Define glioma cell lines and make a separte column for it
depmap_cellinfo = read.csv("CRISPR/DepMap_22Q1/sample_info.csv")
depmap_cellinfo <- depmap_cellinfo[,c('DepMap_ID',"primary_disease","lineage_subtype",'lineage_sub_subtype',
                                      'cell_line_name','primary_or_metastasis','RRID')]
depmap_cellinfo$Glioma <- ''
depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma','Glioma'] <- depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma','lineage_sub_subtype']
#depmap_cellinfo <- depmap_cellinfo[depmap_cellinfo$lineage_subtype =='glioma',] 
#depmap_cellinfo <- depmap_cellinfo[depmap_cellinfo$lineage_sub_subtype!='',] 
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-001198','cell_line_name'] <- 'SNB19'
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-001000','cell_line_name'] <- '1321N1'
depmap_cellinfo[depmap_cellinfo$DepMap_ID == 'ACH-000591','cell_line_name'] <- 'LN235'
depmap_cellinfo$Subtype <- ''
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="astrocytoma" ] <-'AST'
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="glioblastoma" ] <-'GBM'
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_sub_subtype=="oligodendroglioma" ] <-'ODG'

intersect(str_to_lower(vitro$Cell_lines),str_to_lower(depmap_cellinfo$cell_line_name))
intersect(vitro$Cell_lines,depmap_cellinfo$cell_line_name)
vitro$Cell_lines_ <- str_to_lower(vitro$Cell_lines)
depmap_cellinfo$cell_line_name_ <- str_to_lower(depmap_cellinfo$cell_line_name)
depmap_cellinfo$cell_line_name <- str_replace_all(depmap_cellinfo$cell_line_name,'-','')
depmap_cellinfo$cell_line_name <- str_replace_all(depmap_cellinfo$cell_line_name,' ','')
depmap_cellinfo <- distinct(depmap_cellinfo)

non_gliomas <- c("medulloblastoma","ATRT", "neuroblastoma", "Ewing_sarcoma","meningioma","PNET")
depmap_cellinfo$Subtype[depmap_cellinfo$lineage_subtype %in% non_gliomas] <- 'Non-glioma'

depmap_cellinfo$Subtype[depmap_cellinfo$Cell_lines %in% c('U251MG','U251','U138MG','U373','U118MG','U373MG') ] <-'AST'
depmap_cellinfo$Subtype[depmap_cellinfo$Cell_lines %in% c('C6','U87MG','U87','A172','LN18','LN229','SF539','SNB75','SNB78','T98','U343MG') ] <-'GBM'

colnames(depmap_cellinfo)

# Merge the 2 subtype information from cellosaurus
depmap_cellinfo <- left_join(depmap_cellinfo,df4[,c('RRID','Subtype_cellosaurus')])
depmap_cellinfo$Subtype[depmap_cellinfo$Subtype==''] <- depmap_cellinfo$Subtype_cellosaurus[depmap_cellinfo$Subtype=='']
write_csv(depmap_cellinfo,"CRISPR/DepMap_22Q1/sample_info_braincancer.csv")
table(depmap_cellinfo$Subtype)
table(df4$Subtype_cellosaurus)
