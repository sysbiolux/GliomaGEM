%%% Define toxic essential genes, drug and drug combinations using the
%%% control model
solverOK=changeCobraSolver('ibm_cplex','all');
addpath(genpath("./GliomGEM"));

functions = {'biomass_reaction','biomass_maintenance','DM_atp_c_'};

sko_drugs = readtable('./Integrated_Drug_Db/SKO_result_with_Effective_Targets.csv');
sko_drugs = unique(sko_drugs.Drugs);
dko_drugs = readtable('./Integrated_Drug_Db/DKO_result_with_Effective_Targets.csv');
dko_drugs = unique(dko_drugs.Drugs);

model_recon3d = load('./Generic_Models/Recon3DModel_301.mat');
model_recon3d = model_recon3d.Recon3DModel;

%% Loading the ctrl model
model_ctl = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_CTRL.mat');
model_ctl = model_ctl.context_model_w_transcript;

load('./Generic_Models/dico_short.mat');

integrated_db  = readtable('Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv');
integrated_db = integrated_db(:,1:8);
integrated_db.Properties.VariableNames = {'DrugName','uniprot_id','database','ENTREZ','SYMBOL','clinical_phase','moa','disease_area'};
integrated_db.ENTREZ = string(integrated_db.ENTREZ);
integrated_db = integrated_db(:,{'DrugName','ENTREZ'});

comb_db = load("./Integrated_Drug_Db/Approved_Glioma_Comb.mat");
comb_db =comb_db.comb_db;

sko_results = table();
sko_results.Drugs =  sko_drugs;
sko_results.biomass_reaction = repmat(-1,numel(sko_drugs),1);
sko_results.biomass_maintenance = repmat(-1,numel(sko_drugs),1);
sko_results.DM_atp_c_ = repmat(-1,numel(sko_drugs),1);

dko_results = table();
dko_results.Drugs =  dko_drugs;
dko_results.biomass_reaction = repmat(-1,numel(dko_drugs),1);
dko_results.biomass_maintenance = repmat(-1,numel(dko_drugs),1);
dko_results.DM_atp_c_ = repmat(-1,numel(dko_drugs),1);

for i=1:numel(functions)
    model_ctl = changeObjective(model_ctl,functions{i},1);
    %%%%%% Single drug deletion with the predicted single drugs
    [grRatio_ctl, grRateKO, grRateWT, hasEffect, delRxns_ctl, fluxSolution_ctl] = DrugDeletion_v3(...
    model_ctl, 'FBA', sko_drugs,integrated_db);
    %%%%%% Double drug deletion with the predicted single drugs
    [grRatio_ctl_dko, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution_ctl_dko] = DrugDeletion_v3(...
    model_ctl, 'FBA', dko_drugs,comb_db);
    sko_results{:,1+i} = grRatio_ctl;
    dko_results{:,1+i} = grRatio_ctl_dko;
end

writetable(sko_results,'./Integrated_Drug_Db/SKO_result_CTRL_model.csv');
writetable(dko_results,'./Integrated_Drug_Db/DKO_result_CTRL_model.csv');
