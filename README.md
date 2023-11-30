### GliomaGEM: Glioma metabolic modelling for repurposing of single and combination drugs.###

# Data: expression data for TCGA-GBM and TCGA-LGG from three sources
* CRISPR: DepMap 22Q1 CRISPR screen
* Primary_target_databases: Three drug-target databases (DrugBank,  PROMISCUOUS, Drug Repurposing Hub)
* Primary_target_databases/Orphandrug_GlialTumor.csv : Investigational anti-glioma drugs downloaded for OrhanDrug
* Primary_target_databases/Approved_anti_brain_cancers.csv
* Supplementary File 2.xlsx: Literature data for gene KO, drug response and clinical trials

* Sample_models/Subtypes_Models_w_Transcripts: Built sample models for the three glioma subtypes
* Concensus_models : Concensus models for the three glioma subtypes and the control model

## Data preprocessing ##
## Download the complete metadata for TCGA-GBM and TCGA-LGG
TCGA_Biolinks_Download.R

## PCA of the expression data for the three subtypes
Preprocess_TCGA.R

Jaccard_Similarity.m
Merge_metabolic_genes.m

### Model Building ##
## Generate a consistent Human 1 model
Get_Consistent_Human1.m
## Sample model building over 3 reconstructions, medium/ no medium, and 3 datasets
Model_Building_SampleLevel_Three_Models_CSF.m
## Clustering reaction presence for choosing optimal combination of: data, medium and reconstruction
TCGA_Data_Model_Selection.R

## Concensus model building
Model_Building_Subtype_Three_Models_CSF.m

## Copy best model in sample clustering to Concensus_models
mkdir Concensus_models
cp Sample_models/Subtypes_Models_w_Transcripts/Recon3D_Rahman2015__CSF_Thiele2020/* Concensus_models/

Models_Statistics_Consensus.m

## Single gene deletion and FVA as model quality contro
Model_Gene_Deletion_and_FVA.m
## Prepare depmap CRISPR data to compare approved vs predicted drug targets
Classify_Glioma_Celllines.R
CRISPR_DepMap_Data_Splitting.m

### Drug Deletion ##
## Parse DrugBank V5 from .xml to csv
DrugBank-Parsing.ipynb
## Merge drug-target databases prior the drug deletion
Merge_target_databases.R
## Drug deletion of the three glioma models
Drug_Deletion_Glioma.m
## Toxicity analysi of predicted drugs in the control model
Drug_Deletion_Control.m

## Merce drug responses form HTS for IC50, viability reduction, andt xenograft data

Merge_Drug_Response_Databases_IC50.R
Merge_Drug_Response_Databases_with_Viability.R
Merge_Drug_Response_Databases_Xenografts.R

## Create all main and suppl. figures of the paper
Visualizing_the_results.R
