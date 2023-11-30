%% Single gene deletion
solverOK=changeCobraSolver('ibm_cplex','all');
addpath(genpath("./GliomGEM"));

context_model_ast = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_AST_IDH_mut.mat');
context_model_gbm = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_GBM_IDH_wt.mat');
context_model_odg = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_ODG_IDH_mut_Codel.mat');
context_model_ast = context_model_ast.context_model_w_transcript;
context_model_gbm = context_model_gbm.context_model_w_transcript;
context_model_odg = context_model_odg.context_model_w_transcript;
context_model_ast = changeObjective(context_model_ast,'biomass_reaction');
context_model_gbm = changeObjective(context_model_gbm,'biomass_reaction');
context_model_odg = changeObjective(context_model_odg,'biomass_reaction');

threshold = 0.5;
[grRatio, grRateKO, grRateWT, ~, ~, ~, geneList]= singleGeneDeletion_rFASTCORMICS(...
    context_model_ast, 'FBA',[],0,1);
essential_genes_ast =  geneList(grRatio<= threshold);
[grRatio, grRateKO, grRateWT, ~, ~, ~, geneList]= singleGeneDeletion_rFASTCORMICS(...
    context_model_gbm, 'FBA',[],0,1);
essential_genes_gbm =  geneList(grRatio<= threshold);
[grRatio, grRateKO, grRateWT, ~, ~, ~, geneList]= singleGeneDeletion_rFASTCORMICS(...
    context_model_odg, 'FBA',[],0,1);
essential_genes_odg =  geneList(grRatio<= threshold);
save("./Concensus_models/Essentiality_result_.mat");
load("./Concensus_models/Essentiality_result_.mat");
all_genes_union = [essential_genes_ast;essential_genes_gbm;essential_genes_odg]
all_genes_union = unique(string(all_genes_union));
all_genes_union=rmmissing(all_genes_union);
UpSet_table_genes = zeros(numel(["AST";"GBM";"ODG"]),numel(all_genes_union));
UpSet_table_genes(1,ismember(all_genes_union,essential_genes_ast)) = 1;
UpSet_table_genes(2,ismember(all_genes_union,essential_genes_gbm)) = 1;
UpSet_table_genes(3,ismember(all_genes_union,essential_genes_odg)) = 1;
UpSet_table_genes = array2table(UpSet_table_genes);
UpSet_table_genes.Properties.VariableNames = all_genes_union;
UpSet_table_genes.Properties.RowNames =["____AST";"____GBM";"____ODG"];
writetable(UpSet_table_genes,'./Concensus_models/Essential_genes_UpSet_table_.csv','WriteRowNames',true);


%%% FVA on the context models %%%%
model_recon3d = load('./Generic_Models/Recon3DModel_301.mat');
model_recon3d = model_recon3d.Recon3DModel;
model_recon2 = load('./Generic_Models/consistRecon2_4.mat');
model_recon2 = model_recon2.model;
model_recon2.description = "recon2_model";

optional_settings.func = {'biomass_reaction','DM_atp_c_'}; % forced additional reactions into the  model
ex_rxns_ast = findEX_Rxns_fastcormics(context_model_ast,{'biomass_reaction'},optional_settings.func);
[minFlux_ast,maxFlux_ast] = fluxVariability(context_model_ast,100,'rxnNameList',ex_rxns_ast); % performing FVA

ex_rxns_gbm = findEX_Rxns_fastcormics(context_model_gbm,{'biomass_reaction'},optional_settings.func);
[minFlux_gbm,maxFlux_gbm] = fluxVariability(context_model_gbm,100,'rxnNameList',ex_rxns_gbm); % performing FVA

ex_rxns_odg = findEX_Rxns_fastcormics(context_model_odg,{'biomass_reaction'},optional_settings.func);
[minFlux_odg,maxFlux_odg] = fluxVariability(context_model_odg,100,'rxnNameList',ex_rxns_odg); % performing FVA

ex_rxns = ex_rxns_ast;
medium_ex_T_ast = table(ex_rxns,minFlux_ast,maxFlux_ast);
ex_rxns = ex_rxns_gbm;
medium_ex_T_gbm = table(ex_rxns,minFlux_gbm,maxFlux_gbm);
ex_rxns = ex_rxns_odg;
medium_ex_T_odg = table(ex_rxns,minFlux_odg,maxFlux_odg);

medium_ex_T_union = outerjoin(medium_ex_T_ast,medium_ex_T_gbm,'MergeKeys',true)
medium_ex_T_union = outerjoin(medium_ex_T_union,medium_ex_T_odg,'MergeKeys',true)

%Select only ex from the CSF
medium_rxns = findRxnsFromMets(model_recon3d,medium_df_recon);
medium_ex_T_union = medium_ex_T_union(ismember(medium_ex_T_union.ex_rxns,medium_rxns),:);

medium_ex_formulas = printRxnFormula(model_recon3d,'rxnAbbrList',medium_ex_T_union.ex_rxns,'metNameFlag',true)
medium_ex_T_union.Formula = medium_ex_formulas;
writetable(medium_ex_T_union,'./Concensus_models/Subtypes_Models_FVA.csv','WriteRowNames',true);

%%% FVA on the context models for all reactions %%%%
optional_settings.func = {'biomass_reaction','DM_atp_c_'}; % forced additional reactions into the  model
[minFlux_ast,maxFlux_ast] = fluxVariability(context_model_ast,100); % performing FVA
[minFlux_gbm,maxFlux_gbm] = fluxVariability(context_model_gbm,100); % performing FVA
[minFlux_odg,maxFlux_odg] = fluxVariability(context_model_odg,100); % performing FVA

ex_rxns = ex_rxns_ast;
medium_ex_T_ast = table(ex_rxns,minFlux_ast,maxFlux_ast);
ex_rxns = ex_rxns_gbm;
medium_ex_T_gbm = table(ex_rxns,minFlux_gbm,maxFlux_gbm);
ex_rxns = ex_rxns_odg;
medium_ex_T_odg = table(ex_rxns,minFlux_odg,maxFlux_odg);

medium_ex_T_ast = table(context_model_ast.rxns,minFlux_ast,maxFlux_ast);
medium_ex_T_gbm = table(context_model_gbm.rxns,minFlux_gbm,maxFlux_gbm);
medium_ex_T_odg = table(context_model_odg.rxns,minFlux_odg,maxFlux_odg);
medium_ex_T_ast.Properties.VariableNames = {'Rxn','min_AST','max_AST'};
medium_ex_T_gbm.Properties.VariableNames = {'Rxn','min_GBM','max_GBM'};
medium_ex_T_odg.Properties.VariableNames = {'Rxn','min_ODG','max_ODG'};

medium_ex_T_union = outerjoin(medium_ex_T_ast,medium_ex_T_gbm,'MergeKeys',true)
medium_ex_T_union = outerjoin(medium_ex_T_union,medium_ex_T_odg,'MergeKeys',true)

medium_ex_formulas = printRxnFormula(model_recon3d,'rxnAbbrList',medium_ex_T_union.Rxn,'metNameFlag',true);
medium_ex_T_union.Formula = medium_ex_formulas;
writetable(medium_ex_T_union,'./Concensus_models/Subtypes_Models_FVA_All_Reactions.csv','WriteRowNames',true);
