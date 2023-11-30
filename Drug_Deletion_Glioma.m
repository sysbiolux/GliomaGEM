%% Single drug deletion for approved drugs for glioma
cd("./GliomaGEM/");
addpath(genpath("./GliomGEM"));
changeCobraSolver('ibm_cplex')
%glioma_drugs = unique(["carmustine","lomustine","temozolomide","bevacizumab",...
 %   "everolimus","belzutifan"]);
load('./Generic_Models/dico_short.mat');

integrated_db  = readtable('Integrated_Drug_Db/DrugBank_PROMISCUOUS_DRH.csv');
integrated_db = integrated_db(:,1:8);
integrated_db.Properties.VariableNames = {'DrugName','uniprot_id','database','ENTREZ','SYMBOL','clinical_phase','moa','disease_area'};
integrated_db.ENTREZ = string(integrated_db.ENTREZ);

%glioma_db = integrated_db(ismember(integrated_db.DrugName,glioma_drugs),:);
%glioma_db = unique(glioma_db(:,{'DrugName','ENTREZ'}));

glioma_db_  = readtable('Integrated_Drug_Db/Glioma_drugs_targets.csv');
glioma_db_.Properties.VariableNames = {'DrugName','uniprot_id','database','ENTREZ','SYMBOL','clinical_phase','moa','disease_area'};
glioma_db = unique(glioma_db_(:,{'DrugName','ENTREZ'}));
glioma_db.ENTREZ = string(glioma_db.ENTREZ);

glioma_drugs = unique(glioma_db.DrugName);

context_model_ast = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_AST_IDH_mut.mat');
context_model_gbm = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_GBM_IDH_wt.mat');
context_model_odg = load('./Concensus_models/Recon3D_Rahman2015__CSF_Thiele2020_ODG_IDH_mut_Codel.mat');
context_model_ast = context_model_ast.context_model_w_transcript;
context_model_gbm = context_model_gbm.context_model_w_transcript;
context_model_odg = context_model_odg.context_model_w_transcript;
model_recon3d = load('./Generic_Models/Recon3DModel_301.mat');
model_recon3d = model_recon3d.Recon3DModel;

% Number of shared genes with th drug targets
shared_genes_T = cell2table(cell(numel(glioma_drugs),8));
for i=1:numel(glioma_drugs)
    drug = glioma_drugs(i);
    drug
    drug_targets = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
    genes1 = unique(intersect(strtok(context_model_ast.genes,'.'),drug_targets));
    genes2= unique(intersect(strtok(context_model_gbm.genes,'.'),drug_targets));
    genes3= unique(intersect(strtok(context_model_odg.genes,'.'),drug_targets));
    genes4 = unique(intersect(strtok(model_recon3d.genes,'.'),drug_targets));

    shared_genes_T(i,1) = cellstr(string(numel(genes1)));
    shared_genes_T(i,2) = cellstr(string(numel(genes2)));
    shared_genes_T(i,3) = cellstr(string(numel(genes3)));
    shared_genes_T(i,4) = cellstr(string(numel(genes4)));

    shared_genes_T(i,5) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes1),'SYMBOL'})','; '));
    shared_genes_T(i,6) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes2),'SYMBOL'})','; '));
    shared_genes_T(i,7) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes3),'SYMBOL'})','; '));
    shared_genes_T(i,8) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes4),'SYMBOL'})','; '));

end
shared_genes_T.Drugs = glioma_drugs;
writetable(shared_genes_T,'Integrated_Drug_Db/Glioma_Drugs_in_Models.csv');


%%%%%% Single drug deletion WITH glioma drugs

[grRatio_ast, grRateKO, grRateWT, hasEffect, delRxns_ast, fluxSolution_ast] = DrugDeletion_v3(...
context_model_ast, 'FBA', glioma_drugs,glioma_db);

[grRatio_gbm, grRateKO, grRateWT, hasEffect, delRxns_gbm, fluxSolution_gbm] = DrugDeletion_v3(...
context_model_gbm, 'FBA', glioma_drugs,glioma_db);

[grRatio_odg, grRateKO, grRateWT, hasEffect, delRxns_odg, fluxSolution_odg] = DrugDeletion_v3(...
context_model_odg, 'FBA', glioma_drugs,glioma_db);


save("./Integrated_Drug_Db/SKO_GliomaDrugs_AST.mat","grRatio_ast","delRxns_ast","fluxSolution_ast",'-V7');
save("./Integrated_Drug_Db/SKO_GliomaDrugs_GBM.mat","grRatio_gbm","delRxns_gbm","fluxSolution_gbm",'-V7');
save("./Integrated_Drug_Db/SKO_GliomaDrugs_ODG.mat","grRatio_odg","delRxns_odg","fluxSolution_odg",'-V7');
load("./Integrated_Drug_Db/SKO_GliomaDrugs_AST.mat");
load("./Integrated_Drug_Db/SKO_GliomaDrugs_GBM.mat");
load("./Integrated_Drug_Db/SKO_GliomaDrugs_ODG.mat");

grRatio_T = array2table([grRatio_ast,grRatio_gbm,grRatio_odg]);
grRatio_T.Drugs = glioma_drugs;
grRatio_T.SUM = sum(grRatio_T{:,1:3},2);
grRatio_T = sortrows(grRatio_T,'SUM','descend');
barh(grRatio_T{:,1:3})
ylabel('Drugs')
xlabel('grRatio')
x = 1:numel(grRatio_T.Drugs);
yticks([x])
yticklabels(grRatio_T.Drugs)
set(gca,'YTickLabel',grRatio_T.Drugs);
legend({'Astrocytoma','Glioblastoma','Oligodendroglioma'})
saveas(gcf,'Figure/Fig_S_XX_Approved_drug_deletion.png')

%%%%%%%%%%
% 
% % GBM FBA
% model = context_model_gbm;
% solution = optimizeCbModel(model,'max');
% solution.f %151
% 
% % Carmustine targets
% drug = 'temozolomide';
% drug = 'rosuvastatin';
% drug = 'cannabidiol';
% drug = 'eflornithine';
% drug = 'simvastatin';
% drug ='cyclophosphamide';
% drug = 'ribavirin';
% 
% %drug = 'methotrexate';
% %drug = 'resveratrol';
% model = context_model_gbm;
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
% genes= model.genes(ismember(strtok(model.genes,'.'), genes));
% genes
% [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%      model, cellstr(genes))
% solution = optimizeCbModel(model,'max');
% solution.f
% 
% 
% % Drug combination
% model = context_model_gbm;
% 
% solution = optimizeCbModel(model,'max');
% solution.f
% % Carmustine and methotrexate targets
% %drugs = 'methotrexate';
% drugs = ["resveratrol"]
% drugs = ["ezetimibe","cannabidiol"];
% 
% %drugs = ["carmustine","lomustine","temozolomide","bevacizumab",...
% %    "everolimus","belzutifan"];
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drugs),'ENTREZ'});
% genes= model.genes(ismember(strtok(model.genes,'.'), genes));
% [delModel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%      model, cellstr(genes));
% solution = optimizeCbModel(delModel,'max');
% solution.f
% 
% drugs = ["temozolomide"]
% %drugs = ["carmustine","lomustine","temozolomide","bevacizumab",...
% %    "everolimus","belzutifan"];
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drugs),'ENTREZ'});
% genes= model.genes(ismember(strtok(model.genes,'.'), genes));
% [delModel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%      model, cellstr(genes));
% solution = optimizeCbModel(delModel,'max');
% solution.f
% 
% drugs = ["temozolomide","resveratrol"]
% %drugs = ["carmustine","lomustine","temozolomide","bevacizumab",...
% %    "everolimus","belzutifan"];
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drugs),'ENTREZ'});
% genes= model.genes(ismember(strtok(model.genes,'.'), genes));
% [delModel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%      model, cellstr(genes));
% solution = optimizeCbModel(delModel,'max');
% solution.f
% 
% 
% %% PVC combination
% % Drug combination
% model = context_model_odg;
% 
% solution = optimizeCbModel(model,'max');
% solution.f
% % Carmustine and methotrexate targets
% %drugs = 'methotrexate';
% drugs = ["procarbazine","lomustine","vincristine"]
% %drugs = ["carmustine","lomustine","temozolomide","bevacizumab",...
% %    "everolimus","belzutifan"];
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drugs),'ENTREZ'});
% genes= model.genes(ismember(strtok(model.genes,'.'), genes));
% [delModel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%      model, cellstr(genes));
% solution = optimizeCbModel(delModel,'max');
% solution.f


%%% Create a drug pair database of glioma drugs and approved drugs

approved_db = integrated_db(ismember(integrated_db.clinical_phase,'Launched'),:);
approved_db = unique(approved_db(:,{'DrugName','ENTREZ'}));
approved_db.ENTREZ = string(approved_db.ENTREZ);

approved_drugs = unique(string(approved_db.DrugName));
numel(approved_drugs);

%%%%%% Single drug deletion

[grRatio_ast, grRateKO, grRateWT, hasEffect, delRxns_ast, fluxSolution_ast] = DrugDeletion_v3(...
context_model_ast, 'FBA', approved_drugs,approved_db);

[grRatio_gbm, grRateKO, grRateWT, hasEffect, delRxns_gbm, fluxSolution_gbm] = DrugDeletion_v3(...
context_model_gbm, 'FBA', approved_drugs,approved_db);

[grRatio_odg, grRateKO, grRateWT, hasEffect, delRxns_odg, fluxSolution_odg] = DrugDeletion_v3(...
context_model_odg, 'FBA', approved_drugs,approved_db);


save("./Integrated_Drug_Db/SKO_AST.mat","grRatio_ast","delRxns_ast","fluxSolution_ast",'-V7');
save("./Integrated_Drug_Db/SKO_GBM.mat","grRatio_gbm","delRxns_gbm","fluxSolution_gbm",'-V7');
save("./Integrated_Drug_Db/SKO_ODG.mat","grRatio_odg","delRxns_odg","fluxSolution_odg",'-V7');

load("./Integrated_Drug_Db/SKO_AST.mat");
load("./Integrated_Drug_Db/SKO_GBM.mat");
load("./Integrated_Drug_Db/SKO_ODG.mat");

drugs_sing_ast = approved_drugs(grRatio_ast<0.5);
drugs_sing_gbm = approved_drugs(grRatio_gbm<0.5);
drugs_sing_odg = approved_drugs(grRatio_odg<0.5);

grRatio_ast_ = grRatio_ast(grRatio_ast<0.5);
grRatio_gbm_ = grRatio_gbm(grRatio_gbm<0.5);
grRatio_odg_ = grRatio_odg(grRatio_odg<0.5);

delRxns_ast_ = delRxns_ast(grRatio_ast<0.5);
delRxns_gbm_ = delRxns_gbm(grRatio_gbm<0.5);
delRxns_odg_ = delRxns_odg(grRatio_odg<0.5);

delRxns_ast__ = [];
for i=1:numel(delRxns_ast_)
    delRxns_ast__ =[delRxns_ast__;string(strjoin(delRxns_ast_{i},'; '))];
end
delRxns_gbm__ = [];
for i=1:numel(delRxns_gbm_)
    delRxns_gbm__ =[delRxns_gbm__;string(strjoin(delRxns_gbm_{i},'; '))];
end
    
delRxns_odg__ = [];
for i=1:numel(delRxns_odg_)
    delRxns_odg__ =[delRxns_odg__;string(strjoin(delRxns_odg_{i},'; '))];
end

drugs_sing_ast_T = [drugs_sing_ast,repmat('AST',numel(drugs_sing_ast),1),grRatio_ast_,delRxns_ast__];
drugs_sing_gbm_T = [drugs_sing_gbm,repmat('GBM',numel(drugs_sing_gbm),1),grRatio_gbm_,delRxns_gbm__];
drugs_sing_odg_T = [drugs_sing_odg,repmat('ODG',numel(drugs_sing_odg),1),grRatio_odg_,delRxns_odg__];

drugs_sing_all_T = [drugs_sing_ast_T;drugs_sing_gbm_T;drugs_sing_odg_T];
drugs_sing_all_T = array2table(drugs_sing_all_T);
drugs_sing_all_T.Properties.VariableNames = {'Drugs','Subtype','grRatio','DelRxns'};
writetable(drugs_sing_all_T,'Integrated_Drug_Db/SKO_result.csv');

drugs_sing_all_T = readtable('Integrated_Drug_Db/SKO_result.csv');

% Number of shared genes with the drug targets
sko_drugs = string(unique(drugs_sing_all_T{:,1})');
shared_genes_T = cell2table(cell(numel(sko_drugs),8));
for i=1:numel(sko_drugs)
    drug = sko_drugs(i);
    drug;
    drug_targets = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
    genes1 = unique(intersect(strtok(context_model_ast.genes,'.'),drug_targets));
    genes2= unique(intersect(strtok(context_model_gbm.genes,'.'),drug_targets));
    genes3= unique(intersect(strtok(context_model_odg.genes,'.'),drug_targets));
    genes4 = unique(intersect(strtok(model_recon3d.genes,'.'),drug_targets));

    shared_genes_T(i,1) = cellstr(string(numel(genes1)));
    shared_genes_T(i,2) = cellstr(string(numel(genes2)));
    shared_genes_T(i,3) = cellstr(string(numel(genes3)));
    shared_genes_T(i,4) = cellstr(string(numel(genes4)));

    shared_genes_T(i,5) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes1),'SYMBOL'})','; '));
    shared_genes_T(i,6) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes2),'SYMBOL'})','; '));
    shared_genes_T(i,7) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes3),'SYMBOL'})','; '));
    shared_genes_T(i,8) = cellstr(strjoin(unique(dico{ismember(dico.ENTREZ,genes4),'SYMBOL'})','; '));

end
shared_genes_T.Drugs = sko_drugs';
writetable(shared_genes_T,'Integrated_Drug_Db/SKO_Drugs_in_Models.csv');

%%%%%%%% Identify deleted reactions whose single reaction deletion would
%%%%%%%% close the biomass

del_rxns  = drugs_sing_all_T.DelRxns;
del_rxns_list = [];
for i = 1:numel(del_rxns)
    del_rxns_list = [del_rxns_list;split(del_rxns{i},'; ')];
end
del_rxns_list = unique(del_rxns_list);

solution_ast = optimizeCbModel(context_model_ast,'max');
solution_ast.f
solution_gbm= optimizeCbModel(context_model_gbm,'max');
solution_gbm.f
solution_odg = optimizeCbModel(context_model_odg,'max');
solution_odg.f

rxn_del_T = table(del_rxns_list);
rxn_del_T.AST = repmat("NA",numel(del_rxns_list),1);
rxn_del_T.GBM = repmat("NA",numel(del_rxns_list),1);
rxn_del_T.ODG = repmat("NA",numel(del_rxns_list),1);

rxns_idx = find(ismember(del_rxns_list,context_model_ast.rxns));
x = solution_ast.x(ismember(context_model_ast.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'AST'} = string(num2str(x));

rxns_idx = find(ismember(del_rxns_list,context_model_gbm.rxns));
x = solution_gbm.x(ismember(context_model_gbm.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'GBM'} = string(num2str(x));

rxns_idx = find(ismember(del_rxns_list,context_model_odg.rxns));
x = solution_odg.x(ismember(context_model_odg.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'ODG'} = string(num2str(x));

writetable(rxn_del_T,'Integrated_Drug_Db/SKO_Drugs_Del_Rxns.csv');

% for i = 1:numel(del_rxns_list)
% rxn = del_rxns_list{i};
% [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(...
%     context_model_ast,'FBA',string(rxn),0);

%%%%%% Double drug deletion
comb_db = [];
% Exclude drugs identified by SKO from the approved drugs
approved_drugs_ = setdiff(approved_drugs,drugs_sing_all_T.Drugs);
%glioma_drugs= {"cyclophosphamide"};
for i=1:numel(glioma_drugs)
    glioma_drug = glioma_drugs{i};
    glioma_targets = unique(glioma_db{ismember(glioma_db.DrugName,glioma_drug),'ENTREZ'});
    i
    for k=1:numel(approved_drugs_)
        approved_drug = approved_drugs_{k};
        approved_targets = unique(approved_db{ismember(approved_db.DrugName,approved_drug),'ENTREZ'});
        targets_union = union(glioma_targets,approved_targets);
        if size(targets_union,2)>1
            targets_union = targets_union';
        end
       
        comb_name = strcat(glioma_drug,' ; ',approved_drug);
        comb_t = [repmat(comb_name,numel(targets_union),1),targets_union];
        comb_db = [comb_db;comb_t];
    end
end
comb_db = array2table(comb_db);
comb_db.Properties.VariableNames = {'DrugName','ENTREZ'};

save("./Integrated_Drug_Db/Approved_Glioma_Comb.mat","comb_db",'-V7');
comb_db = load("./Integrated_Drug_Db/Approved_Glioma_Comb.mat");
comb_db =comb_db.comb_db;

drugs_comb = unique(comb_db.DrugName);
tic

[grRatio_ast, grRateKO, grRateWT, hasEffect, delRxns_ast, fluxSolution_ast] = DrugDeletion_v3(...
context_model_ast, 'FBA', drugs_comb,comb_db);


[grRatio_gbm, grRateKO, grRateWT, hasEffect, delRxns_gbm, fluxSolution_gbm] = DrugDeletion_v3(...
context_model_gbm, 'FBA', drugs_comb,comb_db);


[grRatio_odg, grRateKO, grRateWT, hasEffect, delRxns_odg, fluxSolution_odg] = DrugDeletion_v3(...
context_model_odg, 'FBA', drugs_comb,comb_db);

toc

save("./Integrated_Drug_Db/DKO_AST.mat","grRatio_ast","delRxns_ast","fluxSolution_ast",'-V7');
save("./Integrated_Drug_Db/DKO_GBM.mat","grRatio_gbm","delRxns_gbm","fluxSolution_gbm",'-V7');
save("./Integrated_Drug_Db/DKO_ODG.mat","grRatio_odg","delRxns_odg","fluxSolution_odg",'-V7');

%,"grRatio_ast","delRxns_ast","fluxSolution_ast"
load("./Integrated_Drug_Db/DKO_AST.mat");
load("./Integrated_Drug_Db/DKO_GBM.mat");
load("./Integrated_Drug_Db/DKO_ODG.mat");

drugs_comb_ast = drugs_comb(grRatio_ast<0.5);
drugs_comb_gbm = drugs_comb(grRatio_gbm<0.5);
drugs_comb_odg = drugs_comb(grRatio_odg<0.5);

grRatio_ast_ = grRatio_ast(grRatio_ast<0.5);
grRatio_gbm_ = grRatio_gbm(grRatio_gbm<0.5);
grRatio_odg_ = grRatio_odg(grRatio_odg<0.5);


delRxns_ast_ = delRxns_ast(grRatio_ast<0.5);
delRxns_gbm_ = delRxns_gbm(grRatio_gbm<0.5);
delRxns_odg_ = delRxns_odg(grRatio_odg<0.5);

delRxns_ast__ = [];
for i=1:numel(delRxns_ast_)
    delRxns_ast__ =[delRxns_ast__;string(strjoin(delRxns_ast_{i},'; '))];
end
delRxns_gbm__ = [];
for i=1:numel(delRxns_gbm_)
    delRxns_gbm__ =[delRxns_gbm__;string(strjoin(delRxns_gbm_{i},'; '))];
end
    
delRxns_odg__ = [];
for i=1:numel(delRxns_odg_)
    delRxns_odg__ =[delRxns_odg__;string(strjoin(delRxns_odg_{i},'; '))];
end

drugs_comb_ast_T = [drugs_comb_ast,repmat('AST',numel(drugs_comb_ast),1),grRatio_ast_,delRxns_ast__];
drugs_comb_gbm_T = [drugs_comb_gbm,repmat('GBM',numel(drugs_comb_gbm),1),grRatio_gbm_,delRxns_gbm__];
drugs_comb_odg_T = [drugs_comb_odg,repmat('ODG',numel(drugs_comb_odg),1),grRatio_odg_,delRxns_odg__];

drugs_comb_all_T = [drugs_comb_ast_T;drugs_comb_gbm_T;drugs_comb_odg_T];
drugs_comb_all_T = array2table(drugs_comb_all_T);
drugs_comb_all_T.Properties.VariableNames = {'Drugs','Subtype','grRatio','DelRxns'};

writetable(drugs_comb_all_T,'Integrated_Drug_Db/DKO_Combination_result.csv');

%%%%%%%% Identify deleted reactions whose single reaction deletion would
%%%%%%%% close the biomass

del_rxns  = drugs_comb_all_T.DelRxns;
del_rxns_list = [];
for i = 1:numel(del_rxns)
    del_rxns_list = [del_rxns_list;split(del_rxns{i},'; ')];
end
del_rxns_list = unique(del_rxns_list);

solution_ast = optimizeCbModel(context_model_ast,'max');
solution_ast.f
solution_gbm= optimizeCbModel(context_model_gbm,'max');
solution_gbm.f
solution_odg = optimizeCbModel(context_model_odg,'max');
solution_odg.f

rxn_del_T = table(del_rxns_list);
rxn_del_T.AST = repmat("NA",numel(del_rxns_list),1);
rxn_del_T.GBM = repmat("NA",numel(del_rxns_list),1);
rxn_del_T.ODG = repmat("NA",numel(del_rxns_list),1);

rxns_idx = find(ismember(del_rxns_list,context_model_ast.rxns));
x = solution_ast.x(ismember(context_model_ast.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'AST'} = string(num2str(x));

rxns_idx = find(ismember(del_rxns_list,context_model_gbm.rxns));
x = solution_gbm.x(ismember(context_model_gbm.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'GBM'} = string(num2str(x));

rxns_idx = find(ismember(del_rxns_list,context_model_odg.rxns));
x = solution_odg.x(ismember(context_model_odg.rxns,del_rxns_list(rxns_idx)));
rxn_del_T{rxns_idx,'ODG'} = string(num2str(x));

writetable(rxn_del_T,'Integrated_Drug_Db/DKO_Drugs_Del_Rxns.csv');


%% The targets of 3 drugs predicetd by DKO
drugs = ["urea","fludarabine","ganciclovir"] ;
drugs_db = integrated_db(ismember(integrated_db.DrugName,drugs),:);


%% Urea + anticancer combination
% Drug combination
model = context_model_gbm;

solution = optimizeCbModel(model,'max');
solution.f
%drugs  = ['ganciclovir'];
%drugs = ['fludarabine'];
drugs = ["fludarabine","bevacizumab"]
%drugs = ["carmustine","lomustine","temozolomide","bevacizumab",...
%    "everolimus","belzutifan"];
genes = unique(integrated_db{ismember(integrated_db.DrugName,drugs),'ENTREZ'});
genes= model.genes(ismember(strtok(model.genes,'.'), genes));
[delModel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
     model, cellstr(genes));
solution = optimizeCbModel(delModel,'max');
solution.f

% % Carmustine targets
% %drug = 'carmustine';
% drug = 'methotrexate';
% 
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
% genes= context_model_odg.genes(ismember(strtok(context_model_odg.genes,'.'), genes));
% 
% solution = optimizeCbModel(context_model_odg,'max');
% solution.f
% 
% [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
%     context_model_odg, cellstr(genes))
% solution = optimizeCbModel(model,'max');
% solution.f
% 
% genes = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
% genes= context_model_odg.genes(ismember(strtok(context_model_odg.genes,'.'), genes));
% 
% [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(...
%     context_model_odg, genes)
% solution = optimizeCbModel(model,'max');
% solution.f
% 
% 
% Single Reaction deletion
[grRatio_r_ast, grRateKO, grRateWT, hasEffect, delRxn_r_ast, fluxSolution_r_ast]=singleRxnDeletion(...
    context_model_ast);
[grRatio_r_gbm, grRateKO, grRateWT, hasEffect, delRxn_r_gbm, fluxSolution_r_gbm]=singleRxnDeletion(...
    context_model_gbm);
[grRatio_r_odg, grRateKO, grRateWT, hasEffect, delRxn_r_odg, fluxSolution_r_odg]=singleRxnDeletion(...
    context_model_odg);

rxn_ast_t = array2table([delRxn_r_ast,string(num2str(grRatio_r_ast))]);
rxn_gbm_t = array2table([delRxn_r_gbm,string(num2str(grRatio_r_gbm))]);
rxn_odg_t = array2table([delRxn_r_odg,string(num2str(grRatio_r_odg))]);
rxn_ast_t.Properties.VariableNames={'Rxn','grRatio_AST'};
rxn_gbm_t.Properties.VariableNames={'Rxn','grRatio_GBM'};
rxn_odg_t.Properties.VariableNames={'Rxn','grRatio_ODG'};

rxn_T = outerjoin(rxn_gbm_t,rxn_ast_t,'Keys',{'Rxn'},'MergeKeys',true);
rxn_T = outerjoin(rxn_T,rxn_odg_t,'Keys',{'Rxn'},'MergeKeys',true);
writetable(rxn_T,'Integrated_Drug_Db/Single_Reaction_Deletion.csv','WriteRowNames',true);


%%% Backword selection of the effective genes from the target list

drugs_sing_all_T = readtable('Integrated_Drug_Db/SKO_result.csv');
drugs_sing_all_T.All_Targets = repmat("0",numel(drugs_sing_all_T.Drugs ),1);
drugs_sing_all_T.Effective_Targets = repmat("0",numel(drugs_sing_all_T.Drugs ),1);

Subtypes= {'AST','GBM','ODG'};
for s=1:numel(Subtypes)   
    Subtype=Subtypes{s};
    if Subtype=='AST'
        model = context_model_ast;
    elseif Subtype=='GBM'
        model = context_model_gbm;
    else
        model = context_model_odg;
    end
    drugs= drugs_sing_all_T{ismember(drugs_sing_all_T.Subtype,Subtype),'Drugs'};
    solution = optimizeCbModel(model,'max');
    biomass_original = solution.f;
    for d=1:numel(drugs)
        drug = drugs{d};
        genes = unique(integrated_db{ismember(integrated_db.DrugName,drug),'ENTREZ'});
        genes= model.genes(ismember(strtok(model.genes,'.'), genes));
        %Calculate the total biomass
        [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
         model, cellstr(genes));
        solution = optimizeCbModel(delmodel,'max');
        biomass_total = solution.f;

        genes_new = genes;

       for n=1:numel(genes)
           gene= genes{n};
           genes_excluded = setdiff(genes,gene);
           %Calculate the new biomass after removing this gene from the
           %list
            [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
             model, cellstr(genes_excluded));
            solution = optimizeCbModel(delmodel,'max');
            biomass_new = solution.f;

           %Calculate the new biomass after using only this gene
            [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
             model, cellstr(gene));
            solution = optimizeCbModel(delmodel,'max');
            biomass_new_2 = solution.f;
            % if absence of gene in deletion equals the total efffect of
            % targets
            if round(biomass_new,14) == round(biomass_total,14) && round(...
                    biomass_new_2,14) == round(biomass_original,14)
             % AND if the single gene deletion equals the biomass of undel
             % model
                    
                genes_new = setdiff(genes_new,gene);
            end
       end
       idx = ismember(drugs_sing_all_T.Subtype,Subtype) & ismember(drugs_sing_all_T.Drugs,drug);
       drugs_sing_all_T{idx,'All_Targets'} = string(strjoin(genes,'; '));
       drugs_sing_all_T{idx,'Effective_Targets'} = string(strjoin(genes_new,'; '));
    end
end
writetable(drugs_sing_all_T,'Integrated_Drug_Db/SKO_result_with_Effective_Targets.csv','WriteRowNames',true);

%%%% the same for DKO
comb_db = load("./Integrated_Drug_Db/Approved_Glioma_Comb.mat");
comb_db =comb_db.comb_db;

drugs_comb_all_T = readtable('Integrated_Drug_Db/DKO_Combination_result.csv');
drugs_comb_all_T.All_Targets = repmat("0",numel(drugs_comb_all_T.Drugs ),1);
drugs_comb_all_T.Effective_Targets = repmat("0",numel(drugs_comb_all_T.Drugs ),1);

Subtypes= {'AST','GBM','ODG'};
for s=1:numel(Subtypes)   
    Subtype=Subtypes{s};
    if Subtype=='AST'
        model = context_model_ast;
    elseif Subtype=='GBM'
        model = context_model_gbm;
    else
        model = context_model_odg;
    end
    drugs= drugs_comb_all_T{ismember(drugs_comb_all_T.Subtype,Subtype),'Drugs'};
    solution = optimizeCbModel(model,'max');
    biomass_original = solution.f;
    for d=1:numel(drugs)
        drug = drugs{d};
        genes = unique(comb_db{ismember(comb_db.DrugName,drug),'ENTREZ'});
        genes= model.genes(ismember(strtok(model.genes,'.'), genes));
        %Calculate the total biomass
        [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
         model, cellstr(genes));
        solution = optimizeCbModel(delmodel,'max');
        biomass_total = solution.f;

        genes_new = genes;

       for n=1:numel(genes)
           gene= genes{n};
           genes_excluded = setdiff(genes,gene);
           %Calculate the new biomass after removing this gene from the
           %list
            [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
             model, cellstr(genes_excluded));
            solution = optimizeCbModel(delmodel,'max');
            biomass_new = solution.f;

           %Calculate the new biomass after using only this gene
            [delmodel, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes_rFASTCORMICS(...
             model, cellstr(gene));
            solution = optimizeCbModel(delmodel,'max');
            biomass_new_2 = solution.f;
            % if absence of gene in deletion equals the total efffect of
            % targets
            if round(biomass_new,3) == round(biomass_total,3) && round(...
                    biomass_new_2,3) == round(biomass_original,3)
             % AND if the single gene deletion equals the biomass of undel
             % model
                    
                genes_new = setdiff(genes_new,gene);
            end
       end
       idx = ismember(drugs_comb_all_T.Subtype,Subtype) & ismember(drugs_comb_all_T.Drugs,drug);
       drugs_comb_all_T{idx,'All_Targets'} = string(strjoin(genes,'; '));
       drugs_comb_all_T{idx,'Effective_Targets'} = string(strjoin(genes_new,'; '));
    end
end
writetable(drugs_comb_all_T,'Integrated_Drug_Db/DKO_result_with_Effective_Targets.csv','WriteRowNames',true);



%%%%%% Single drug deletion WITH 19 dug predicted in combination
drugs_comb_all_T = readtable('Integrated_Drug_Db/DKO_Combination_result.csv');
drugs_comb = drugs_comb_all_T.Drugs;
drugs_comb = split(drugs_comb,' ;');
drugs_comb = unique([drugs_comb(:,1);drugs_comb(:,2)]);
drugs_comb_db = integrated_db(find(ismember(integrated_db.DrugName,drugs_comb)),:);
[grRatio_ast, grRateKO, grRateWT, hasEffect, delRxns_ast, fluxSolution_ast] = DrugDeletion_v3(...
context_model_ast, 'FBA', drugs_comb,drugs_comb_db);

[grRatio_gbm, grRateKO, grRateWT, hasEffect, delRxns_gbm, fluxSolution_gbm] = DrugDeletion_v3(...
context_model_gbm, 'FBA', drugs_comb,drugs_comb_db);

[grRatio_odg, grRateKO, grRateWT, hasEffect, delRxns_odg, fluxSolution_odg] = DrugDeletion_v3(...
context_model_odg, 'FBA', drugs_comb,drugs_comb_db);
table2array(cell2table(drugs_comb))
drugs_comb_ast_T = [table2array(cell2table(drugs_comb)),string(repmat('AST',numel(drugs_comb),1)),grRatio_ast];
drugs_comb_gbm_T = [table2array(cell2table(drugs_comb)),string(repmat('GBM',numel(drugs_comb),1)),grRatio_gbm];
drugs_comb_odg_T = [table2array(cell2table(drugs_comb)),string(repmat('ODG',numel(drugs_comb),1)),grRatio_odg];

drugs_comb_all_T = [drugs_comb_ast_T;drugs_comb_gbm_T;drugs_comb_odg_T];
drugs_comb_all_T = array2table(drugs_comb_all_T);
drugs_comb_all_T.Properties.VariableNames = {'Drugs','Subtype','grRatio'};

writetable(drugs_comb_all_T,'Integrated_Drug_Db/DKO_Combination_individual_grRatio_result.csv');

%  DKO withouth cut-off on grRatio
drugs_comb_ast = drugs_comb(grRatio_ast<=1);
drugs_comb_gbm = drugs_comb(grRatio_gbm<=1);
drugs_comb_odg = drugs_comb(grRatio_odg<=1);

grRatio_ast_ = grRatio_ast(grRatio_ast<=1);
grRatio_gbm_ = grRatio_gbm(grRatio_gbm<=1);
grRatio_odg_ = grRatio_odg(grRatio_odg<=1);


delRxns_ast_ = delRxns_ast(grRatio_ast<=1);
delRxns_gbm_ = delRxns_gbm(grRatio_gbm<=1);
delRxns_odg_ = delRxns_odg(grRatio_odg<=1);

delRxns_ast__ = [];
for i=1:numel(delRxns_ast_)
    delRxns_ast__ =[delRxns_ast__;string(strjoin(delRxns_ast_{i},'; '))];
end
delRxns_gbm__ = [];
for i=1:numel(delRxns_gbm_)
    delRxns_gbm__ =[delRxns_gbm__;string(strjoin(delRxns_gbm_{i},'; '))];
end
    
delRxns_odg__ = [];
for i=1:numel(delRxns_odg_)
    delRxns_odg__ =[delRxns_odg__;string(strjoin(delRxns_odg_{i},'; '))];
end

drugs_comb_ast_T = [drugs_comb_ast,repmat('AST',numel(drugs_comb_ast),1),grRatio_ast_,delRxns_ast__];
drugs_comb_gbm_T = [drugs_comb_gbm,repmat('GBM',numel(drugs_comb_gbm),1),grRatio_gbm_,delRxns_gbm__];
drugs_comb_odg_T = [drugs_comb_odg,repmat('ODG',numel(drugs_comb_odg),1),grRatio_odg_,delRxns_odg__];

drugs_comb_all_T = [drugs_comb_ast_T;drugs_comb_gbm_T;drugs_comb_odg_T];
drugs_comb_all_T = array2table(drugs_comb_all_T);
drugs_comb_all_T.Properties.VariableNames = {'Drugs','Subtype','grRatio','DelRxns'};

writetable(drugs_comb_all_T,'Integrated_Drug_Db/DKO_Combination_result_all_predictions.csv');
