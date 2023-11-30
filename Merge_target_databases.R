### Merging 3 drug-target databases
### DrugBank , PROMISCOUS, Drug Repurpusing Hub
### Drug approval status was extracted from Drug Repurpusing Hub

##### Modified from https://github.com/sysbiolux/Herbal_drug_prediction/blob/main/scripts/Merging_target_database/Merging_target_databases_Cancer_Drugs.R
##### ‘Identification of Bruceine D as a Drug Candidate against Breast Cancer through Genome-Scale Metabolic Modelling and Cell Viability Assay’ Claudia Cipriani, Maria Pires Pacheco, Ali Kishk, Maryem Wachich, Daniel Abankwa, Elisabeth Schaffner-Reckinger, Thomas Sauter Pharmaceuticals 2022, 15(2), 179 University of Luxembourg.

setwd('./GliomGEM')

library(tidyverse)
library(sjmisc)
library(stringr)

drug_bank <- read_csv("./Primary_target_databases/DrugBank_2022/DrugBank_Drug_Target.csv")

promiscuous_drug_information <- read_delim("./Primary_target_databases/Promiscuous/drug_information.csv",delim='|')
promiscuous_drug_target_interactions <- read_delim("./Primary_target_databases/Promiscuous/drug_target_interactions.csv",delim=';', escape_double = FALSE, trim_ws = TRUE)
colnames(promiscuous_drug_target_interactions)[2:3] <- c('pid','uniprot_id')

colnames(drug_bank)
colnames(promiscuous_drug_information)
colnames(promiscuous_drug_target_interactions)

promiscuous_database <- left_join(promiscuous_drug_target_interactions,promiscuous_drug_information)
promiscuous_database <- distinct(promiscuous_database[,c("name","uniprot_id")])
promiscuous_database$database <- 'promiscuous'

drugbank_database <- distinct(drug_bank[,c('name','uniprot_id')])
drugbank_database$database <- 'DrugBank'


# merge the 2 databases
all_databases <- rbind(drugbank_database,promiscuous_database)

##Map using in house dico disctionary
dico = read_csv('./dico_201911.csv')
dico <- dico[,c('ENTREZ','UniProt','SYMBOL')]
colnames(dico) <- c('entrez_id','uniprot_id','SYMBOL')
dico <- na.omit(dico)

all_databases <- left_join(all_databases,dico)
all_databases <- distinct(na.omit(all_databases))
colnames(all_databases)

repurpus_hub <- read_delim("./Primary_target_databases/repurposing_drugs_20200324.txt",skip = 9,delim = '\t')
colnames(repurpus_hub) <- c('name',"clinical_phase" ,"moa","SYMBOL","disease_area", "indication")
repurpus_hub_database  <- repurpus_hub[,c('name','SYMBOL')]
repurpus_hub_database$database <- "RepurpusingHub"
repurpus_hub_database %>%  separate_rows(SYMBOL,sep  = "\\|") ->repurpus_hub_database
repurpus_hub_database <- left_join(repurpus_hub_database,dico)
repurpus_hub_database <- repurpus_hub_database[,colnames(all_databases)]

# ## Add manually DrugBank 5 glioma drugs
# drugbank_glioma  <- read_delim('./Primary_target_databases/DrugBank_V5_Glioma_Druugs.csv')
# drugbank_glioma <- separate_rows(drugbank_glioma,SYMBOL,sep="; ")
# drugbank_glioma <- drugbank_glioma[,c('name','SYMBOL')]
# drugbank_glioma <- na.omit(left_join(drugbank_glioma,dico))
# drugbank_glioma$name <- str_to_lower(drugbank_glioma$name)
# drugbank_glioma <- left_join(drugbank_glioma,repurpus_hub_approval)
# drugbank_glioma$database <- "DrugBank2022"
# drugbank_glioma <- drugbank_glioma[,colnames(all_databases)]
# all_databases <- rbind(all_databases,drugbank_glioma)


### Merge DRH
all_databases <- rbind(all_databases,repurpus_hub_database)
all_databases$name <- str_to_lower(all_databases$name)
## Extract approval status from Drug Repurpusing Hub
repurpus_hub_approval <- repurpus_hub[,c('name',"clinical_phase" ,"moa","disease_area")]
all_databases <- left_join(all_databases,repurpus_hub_approval)
#Standardize the drug names for investagional antigliomas
all_databases$name <- str_replace(all_databases$name ,"5-fluorouracil","fluorouracil")
all_databases$name <- str_replace(all_databases$name ,"l-eflornithine","eflornithine")

table(all_databases[,c('database')])

write_csv(all_databases, "./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv")
write_csv(all_databases,gzfile( "./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv.gz"))


# ## Read TTD 
# ttd_drugs <-vroom('Primary_target_databases/TTD/P1-03-TTD_crossmatching.txt',skip = 28)
# ttd_drugs %>%
#   pivot_wider(names_from = TTDDRUID,
#               values_from = D00AAN...3) -> ttd_drugs_longer
# ttd_drugs_longer <- ttd_drugs_longer[,c('D00AAN...1','DRUGNAME')]
# colnames(ttd_drugs_longer) <- c('DRUGID','name')
# 
# ttd_targets <-vroom('Primary_target_databases/TTD/P2-01-TTD_uniprot_all.txt',skip = 22)
# ttd_targets %>%
#   pivot_wider(names_from = TARGETID,
#               values_from = T00032...3) -> ttd_targets_longer
# ttd_targets_longer <- ttd_targets_longer[,c('T00032...1','UNIPROID')]
# colnames(ttd_targets_longer) <- c('TARGETID','SYMBOL')
# #ttd_targets_longer$SYMBOL <- str_split(ttd_targets_longer$SYMBOL,'\\(', simplify=T)[,2]
# ttd_targets_longer <- separate_rows(ttd_targets_longer,SYMBOL,sep = '-')
# ttd_targets_longer <- separate_rows(ttd_targets_longer,SYMBOL,sep = '; ')
# ttd_targets_longer <- separate_rows(ttd_targets_longer,SYMBOL,sep = '/')
# ttd_targets_longer  <- ttd_targets_longer[str_detect(ttd_targets_longer$SYMBOL,'HUMAN'),]
# 
# ttd_targets_longer$SYMBOL  <- str_split(ttd_targets_longer$SYMBOL,'\\_', simplify=T)[,1]
# 
# ttd_interact <-vroom('Primary_target_databases/TTD/P1-09-Target_compound_activity.txt')
# colnames(ttd_interact) <- c('TARGETID','DRUGID','Pubchem','Activity')
# 
# ttd_database <- left_join(ttd_interact,ttd_targets_longer)
# ttd_database <- na.omit(left_join(ttd_database,ttd_drugs_longer))
# ttd_database <- na.omit(left_join(ttd_database,dico))
# ttd_database$database <- "TTD"
# ttd_database <- ttd_database[,colnames(all_databases)]
# 
# ### Merge TTD
# all_databases <- rbind(all_databases,ttd_database)
# all_databases$name <- str_to_lower(all_databases$name)
# ## Extract approval status from Drug Repurpusing Hub
# repurpus_hub_approval <- repurpus_hub[,c('name',"clinical_phase" ,"moa","disease_area")]
# all_databases <- left_join(all_databases,repurpus_hub_approval)
# 
# table(all_databases[,c('database')])
# 
# write_csv(all_databases, "./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH_TTD.csv")
# write_csv(all_databases,gzfile( "./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH_TTD.csv.gz"))

## Extract approval drugs from drug bank
#drukbank_approval  <- distinct(drug_bank[,c('name','groups')])
#drukbank_approval <- drukbank_approval[str_detect(drukbank_approval$groups,'^approved'),]
#unique(drukbank_approval$groups)

#all_databases <- left_join(all_databases,drukbank_approval)
#all_databases <- distinct(na.omit(all_databases))

# The approval status for some drugs are not accurate such as Temozolamide

# Find the cancer drugs in the 3 databases that CONTAINS approved and in clinical trials glioma drug names
all_databases<- read_csv("./Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv")
glioma_drugs <- read_csv("./Primary_target_databases/Approved_anti_brain_cancers.csv")
#glioma_drugs <- c("carmustine","lomustine","temozolomide","bevacizumab","everolimus"
#                  ,"belzutifan","doxorubicin","cyclophosphamide","vincr")
glioma_drugs_2 <- read_delim('Primary_target_databases/Orphandrug_GlialTumor.csv',delim='\t')
glioma_drugs_2$Drugs <- str_to_lower(glioma_drugs_2$Drugs)
glioma_drugs_2

queries <- all_databases %>% 
  filter(str_detect(name, paste(glioma_drugs$Drugs, collapse = "|")))
unique(queries$name)
queries <- all_databases %>% 
  filter(str_detect(name, paste(glioma_drugs_2$Drugs, collapse = "|")))
unique(queries$name)

mydrugs_in_all_databases <- distinct(all_databases[all_databases$name %in% glioma_drugs$Drugs,])
mydrugs_in_all_databases_2 <- distinct(all_databases[all_databases$name %in% glioma_drugs_2$Drugs,])

#mydrugs_in_all_databases <- filter(all_databases, str_match(paste(glioma_drugs, collapse="|"), name))
#mydrugs_in_all_databases_2 <- filter(all_databases, str_match(paste(gliosma_drugs_2, collapse="|"), name))

mydrugs_in_all_databases<- rbind(mydrugs_in_all_databases,mydrugs_in_all_databases_2)

mydrugs_in_all_databases <- distinct(mydrugs_in_all_databases) 
mydrugs_in_all_databases <- mydrugs_in_all_databases[!str_detect(mydrugs_in_all_databases$name," \\+ "),]

mydrugs_in_all_databases %>% drop_na(uniprot_id) -> mydrugs_in_all_databases

# Add manually the 2 approved combinations (PCV) and (trametinib-dabrafenib)
comb_brafi <- mydrugs_in_all_databases[mydrugs_in_all_databases$name %in% c("trametinib","dabrafenib"),]
comb_pcv <- mydrugs_in_all_databases[mydrugs_in_all_databases$name %in% c("lomustine","procarbazine","vincristine"),]
comb_brafi$name <- "trametinib-dabrafenib"
comb_pcv$name <- "PCV combination"
mydrugs_in_all_databases <- rbind(mydrugs_in_all_databases,comb_brafi,comb_pcv)

# Number of drugs in each database
mydrugs_in_all_databases %>% distinct(database,name)  %>% count(database)
mydrugs_in_all_databases %>% distinct(name,SYMBOL) %>% count(name)

x <- table(mydrugs_in_all_databases[,c('name','database')])
x
x <- as.data.frame(x)
mydrugs_in_all_databases %>% distinct(database,name)  %>% count(database) -> x$total_count
# The number of drug combined in the 2 databases
mydrugs_in_all_databases %>% distinct(name)  %>% count()

write_csv(mydrugs_in_all_databases, "./Integrated_Drug_Db/Glioma_drugs_targets.csv")
