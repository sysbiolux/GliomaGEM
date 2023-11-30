### Integrating PRISM secondary, GDSC1000, GDSC2000, CTRPv2 for IC50 only
## Drug, cell line, IC50

### PRISM secondary ##
prism_sc <- read_csv('./Drug_Response_Databases/secondary-screen-dose-response-curve-parameters.csv')
colnames(prism_sc)[2] <- 'DepMap_ID'
prism_sc <- prism_sc[prism_sc$passed_str_profiling==TRUE,]
prism_sc <- as.data.frame(prism_sc)

prism_sc <- prism_sc[,c('name','ic50','DepMap_ID','ccle_name')]
#
depmap_cellinfo = read.csv("./Drug_Response_Databases/sample_info.csv")
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

# CTRPv2
ctrp <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt',delim='\t')
colnames(ctrp)
ctrp <- ctrp[,c("master_cpd_id","experiment_id" ,"apparent_ec50_umol")] # area_under_curve
ctrp_cpd <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt',delim='\t')
ctrp_cll <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt',delim='\t')
ctrp_exp <- read_delim('./Drug_Response_Databases/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt',delim='\t')
ctrp_cpd <- ctrp_cpd[,c('master_cpd_id','cpd_name')]
ctrp_cll <- ctrp_cll[,c('master_ccl_id','ccl_name')]
ctrp_exp <- ctrp_exp[,c('experiment_id','master_ccl_id')]
ctrp_cll <- left_join(ctrp_cll,ctrp_exp) # Joining, by = "master_ccl_id"

ctrp <- left_join(ctrp,ctrp_cpd) # Joining, by = "master_cpd_id"
ctrp <- left_join(ctrp,ctrp_cll) # Joining, by = "experiment_id"
ctrp <- ctrp[,c('cpd_name','ccl_name','apparent_ec50_umol')]

### Genentech Cell Line Screening Initiative
gcsi <- read_delim('./Drug_Response_Databases/data_5_Genentech_Cell_Line_Screening_Initiative_(gCSI).csv',delim = ',')
gcsi <- gcsi[,c('Perturbagen','Cell_Line','IC50')]

### Merge the 4 databases
colnames(prism_sc) <- c('Drugs','cell_line_name','IC50')
colnames(gdsc1000) <- c('Drugs','cell_line_name','IC50')
colnames(gdsc2000) <- c('Drugs','cell_line_name','IC50')
colnames(gcsi) <- c('Drugs','cell_line_name','IC50')
colnames(ctrp) <- c('Drugs','cell_line_name','IC50')

prism_sc$Database <-'PRISM Secondary'
gdsc1000$Database <-'GDSC1000'
gdsc2000$Database <-'GDSC2000'
gcsi$Database <-'gCSI'
ctrp$Database <- "CTRPv2"


Merged_ic50 <- rbind(prism_sc,gdsc1000,gdsc2000,gcsi,ctrp)
Merged_ic50 <- distinct(Merged_ic50)
Merged_ic50 <- Merged_ic50[!is.na(Merged_ic50$IC50),]
Merged_ic50 <- Merged_ic50[Merged_ic50$IC50!='Inf',]

Merged_ic50_drugs <-unique(Merged_ic50[,c('Database','Drugs')])

glioma_drugs <- read_csv("./Primary_target_databases/Approved_anti_brain_cancers.csv")
queries <- Merged_ic50_drugs %>% 
  filter(str_detect(Drugs, paste(glioma_drugs$Drugs, collapse = "|")))
unique(queries$Drugs)

write_csv(Merged_ic50,'./Merged_IC50_databases.csv')