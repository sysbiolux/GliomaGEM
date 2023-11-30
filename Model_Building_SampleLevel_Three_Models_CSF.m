HPC =0 ; % paramater for running parfor on HPC
datas = ["Biolinks","RSEM","Rahman2015"];
models = ["Recon2","Recon3D","Human1"];
curation = '_CSF_Thiele2020_DrugMetabolismRemoved';

% avoid building model for recon 2 with _DrugMetabolismRemoved since it
% doesn't have this pathway
if contains(curation,"_DrugMetabolismRemoved")==1
    models = ["Recon3D","Human1"];
end
%curation = 'NoMedium';

%
changeCobraSolver('ibm_cplex');
solverOK=changeCobraSolver('ibm_cplex','all');


load('./Generic_Models/dico_short.mat')
dico_recon3d_ensemble = cellstr(table2cell(dico(:,{'ENSG','ENTREZ'})));
dico_recon3d_genename = cellstr(table2cell(dico(:,{'SYMBOL','ENTREZ'}))); 

%% Load the expression data
tcga_gbm = readtable('./data/TCGABiolinks/TCGA_GBM_FPKM_TCGBiolinks.csv');
tcga_lgg = readtable('./data/TCGABiolinks/TCGA_LGG_FPKM_TCGBiolinks.csv');
intersect(tcga_gbm.Properties.VariableNames,tcga_lgg.Properties.VariableNames);
shared_genes = intersect(string(table2cell(tcga_gbm(:,1))),string(table2cell(tcga_lgg(:,1))));
fpkm = innerjoin(tcga_gbm,tcga_lgg);
rownames_biolinks = cellstr(table2cell(fpkm(:,1)));
fpkm=fpkm(:,2:end);
colnames_biolinks= fpkm.Properties.VariableNames;
fpkm_biolinks = table2array(fpkm);

fpkm = readtable('./data/TCGA_Rahman2015/TCGA_LGG_GBM.csv');
rownames_rahman = cellstr(table2cell(fpkm(:,1)));
fpkm=fpkm(:,2:end);
colnames_rahman= fpkm.Properties.VariableNames;
fpkm_rahman = table2array(fpkm);

fpkm = readtable('./data/TCGA_Ceccarelli2016/TCGA_FPKM.csv');
rownames_rsem = string(table2cell(fpkm(:,1)));
fpkm=fpkm(:,2:end);
colnames_rsem= fpkm.Properties.VariableNames;
fpkm_rsem = table2array(fpkm);

%% Load the consistent generic models
if contains(curation,"_DrugMetabolismRemoved")==1
    model_recon3d = load('./Generic_Models/Recon3DModel_301_NoDrugMetabolism.mat');
    model_recon3d = model_recon3d.model_recon3d_c;
    model_human1 = load('./Generic_Models/Human-GEM_1.5_Consistent_NoDrugMetabolism.mat');
    model_human1 = model_human1.model_human1_c;  
else
    model_recon3d = load('./Generic_Models/Recon3DModel_301.mat');
    model_recon3d = model_recon3d.Recon3DModel;
    model_human1 = load('./Generic_Models/Human-GEM_1.5_Consistent.mat');
    model_human1 = model_human1.ihuman_consistent;
end

model_recon2 = load('./Generic_Models/consistRecon2_4.mat');
model_recon2 = model_recon2.model;
model_recon2.description = "recon2_model";

for i=1:numel(model_human1.rxns)
    model_human1.subSystems{i}=char(model_human1.subSystems{i});
end
model_human1.subSystems = convertStringsToChars(model_human1.subSystems);

curation
numel(model_recon2.rxns)
numel(model_recon3d.rxns)
numel(model_human1.rxns)

%Removing gene version from the model
model_recon3d.genes = string(regexprep(model_recon3d.genes,'\.[0-9]+$',''));
model_recon2.genes = string(regexprep(model_recon2.genes,'\.[0-9]+$',''));

%setting model reconstruction parameters
already_mapped_tag = 0;
consensus_proportion = 0.9;
epsilon = 1e-4;
optional_settings.func = {'DM_atp_c_','biomass_maintenance'};

% inhouse dictionary for recon model
% for other models and data, the user has to create a dictionary using for~
% instance biomart or db2db
dico_recon3d_ensemble = dico(:,{'ENSG','ENTREZ'});
dico_recon3d_genename = dico(:,{'SYMBOL','ENTREZ'});


%% Mapping the CSF medium from Thiele2020 to the 3 models
medium_df = readtable('./Thiele2020_CSF_.csv');
medium_df = unique(medium_df.VMHID);
medium_df_recon = strcat(medium_df,'[e]');
medium_recon2 = intersect(medium_df_recon,model_recon2.mets);
medium_recon3d = intersect(medium_df_recon,model_recon3d.mets);

dico_human1_met = readtable('./Generic_Models/metabolites_09_2021.tsv','FileType','text');
idxs = find(ismember(dico_human1_met.metRecon3DID,medium_df));
medium_human1 = dico_human1_met{idxs,"mets"};
medium_human1 = medium_human1(contains(medium_human1,"s"));
medium_human1 = intersect(medium_human1,model_human1.mets);

%setting model reconstruction parameters
already_mapped_tag = 0;
consensus_proportion = 0.9;

%Recon2
biomass_rxn_recon2 = {'biomass_reaction'};
optional_settings_recon2.func = {'biomass_reaction','DM_atp_c_'}; % forced additional reactions into the  model
%Recon3d
biomass_rxn_recon3d = {'biomass_reaction'};
optional_settings_recon3d.func = {'biomass_reaction','DM_atp_c_'}; % forced additional reactions into the  model
% %Human1 
% optional_settings_human1.medium = medium_human1;
% unpenalizedSystems = {'Transport reactions'};
% unpenalized = model_human1.rxns(ismember(model_human1.subSystems,{'Transport reactions'}));
% optional_settings_human1.unpenalized = unpenalized;
% optional_settings_human1.func = {'biomass_human','EX_atp[e]'};
% not_medium_constrained = 'MAR00696';
% optional_settings_human1.not_medium_constrained = not_medium_constrained;
% biomass_rxn_human1 = {'biomass_human'};

%Human1 
biomass_rxn_human1 = "Generic human cell biomass reaction";%{'biomass_human'};
biomass_rxn_human1 = char(model_human1.rxns(find(ismember(model_human1.rxnNames,biomass_rxn_human1))));
optional_settings_human1.func = {biomass_rxn_human1,'MAR03964'}; %MAR03964 is ATP demand reaction
%%MAR13082

%% Optional setting for medium-constraining
if contains(curation,"_CSF_Thiele2020" ) ==1
    %Recon2
    optional_settings_recon2.medium = medium_recon2;
    unpenalizedSystems = {'Transport, endoplasmic reticular';
        'Transport, extracellular';
        'Transport, golgi apparatus';
        'Transport, mitochondrial';
        'Transport, peroxisomal';
        'Transport, lysosomal';
        'Transport, nuclear'};
    unpenalized = model_recon2.rxns(ismember(model_recon2.subSystems,unpenalizedSystems));
    optional_settings_recon2.unpenalized = unpenalized;
    not_medium_constrained = 'EX_tag_hs(e)';
    optional_settings_recon2.not_medium_constrained = not_medium_constrained;

    %Recon3D
    optional_settings_recon3d.medium = medium_recon3d;
    unpenalizedSystems = {'Transport, endoplasmic reticular';
        'Transport, extracellular';
        'Transport, golgi apparatus';
        'Transport, mitochondrial';
        'Transport, peroxisomal';
        'Transport, lysosomal';
        'Transport, nuclear'};
    unpenalized = model_recon3d.rxns(ismember(string(model_recon3d.subSystems),unpenalizedSystems));
    optional_settings_recon3d.unpenalized = unpenalized;
    not_medium_constrained = 'EX_tag_hs(e)';
    optional_settings_recon3d.not_medium_constrained = not_medium_constrained;
    %Human1
    %Human1
    %optional_settings_human1.medium = medium_human1;
    unpenalizedSystems = {'Transport reactions'};%;'Exchange/demand reactions'};
    x = convertCharsToStrings(model_human1.subSystems);
    %x = model_human1.subSystems;
    unpenalized = model_human1.rxns(find(ismember(x,unpenalizedSystems)));
    optional_settings_human1.unpenalized = unpenalized;
    %not_medium_constrained = {'MAR09153';'MAR09167';'MAR09269';'MAR09154';'MAR10442';'MAR09276';'MAR10492';'MAR09404';'MAR11895'};
    %not_medium_constrained = {'MAR09147';'MAR09153';'MAR09167';'MAR09269';'MAR09276';'MAR09404';'MAR10492';'MAR11895';'MAR13067'};
    % Load not_medium_constarined
    if contains(curation,"_DrugMetabolismRemoved")==1
        not_medium_constrained = load("medium_data/Rxns_zero_flux_Human1_DrugMetabolismRemoved.mat");
    else
        not_medium_constrained = load("medium_data/Rxns_zero_flux_Human1.mat");
    end
    not_medium_constrained = not_medium_constrained.a;
    EX_met_to_open = findMetsFromRxns(model_human1,cellstr(not_medium_constrained))
    medium_human1_extended =[medium_human1;EX_met_to_open];
    optional_settings_human1.medium = medium_human1_extended;
    %optional_settings_human1.not_medium_constrained = '';
end
%% Discretization
discretized_biolinks = discretize_FPKM_skewed(fpkm_biolinks,colnames_biolinks);
%discretized_rsem = discretize_FPKM(fpkm_rsem,colnames_rsem);
%discretized_rahman = discretize_FPKM(fpkm_rahman,colnames_rahman);

% Human1 dico

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 6);
% % % Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["genes", "geneENSTID", "geneENSPID", "geneUniProtID", "geneNames", "geneEntrezID"];
opts.VariableTypes = ["string", "string", "string", "double", "string", "double"];
opts = setvaropts(opts, [2, 3, 5], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 3, 5], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
dico_human = readtable("./Generic_Models/genes_08_2021.tsv", opts);
%ensemble_genes = unique(strsplit(string(dico.geneENSPID),';'));
dico_ens_ens = dico_human(:,{'genes','genes'});
dico_ens_ens =unique(dico_ens_ens,'rows');
dico_ens_ens= table2array(dico_ens_ens);
dico_genename_ens = dico_human(:,{'geneNames','genes'});
dico_genename_ens =unique(dico_genename_ens,'rows');
dico_genename_ens= table2array(dico_genename_ens);
dico_entrez_ens = dico_human(:,{'geneEntrezID','genes'});
dico_entrez_ens =unique(dico_entrez_ens,'rows');
dico_entrez_ens= table2array(dico_entrez_ens);

dico_recon3d_ensemble = dico(:,{'ENSG','ENTREZ'});
dico_recon3d_genename = dico(:,{'SYMBOL','ENTREZ'});
dico_recon3d_genename=table2array(dico_recon3d_genename);
dico_recon3d_ensemble=table2array(dico_recon3d_ensemble);

for d=1:1%numel(datas)
    for m=1:numel(models) 
        %clear generic_name generic_model data optional_settings biomass_rxn dico_ discretized
        generic_name = models{m};
        data = datas{d};
        models_dir = strcat('models/',models{m}, '_TCGA_FPKM_',datas{d},curation,'/');
        if ~exist(models_dir, 'dir')
           mkdir(models_dir)
        end
        
        %% Model parameters
        if generic_name=="Recon2"
           generic_model =  model_recon2;
           optional_settings = optional_settings_recon2;
           biomass_rxn= biomass_rxn_recon2;
        elseif generic_name=="Recon3D"
           generic_model =  model_recon3d;
           optional_settings = optional_settings_recon3d;
           biomass_rxn= biomass_rxn_recon3d;
        elseif generic_name=="Human1"
           generic_model =  model_human1;
           optional_settings = optional_settings_human1;
           biomass_rxn= biomass_rxn_human1;
        end   
        %% Data paarameters
        if data == "Biolinks"
           discretized = discretized_biolinks;
           rownames = rownames_biolinks;
           colnames  = colnames_biolinks;          
        elseif data== "RSEM"
           discretized = discretized_rsem;
           rownames = rownames_rsem;
           colnames  = colnames_rsem;              
        elseif data== "Rahman2015"
           discretized = discretized_rahman;
           rownames = rownames_rahman;
           colnames  = colnames_rahman;  
        end
        
        %% If satament for the dictionary
        if generic_name == "Human1"
            if data =="Biolinks"
                dico_ = dico_ens_ens;
            else
                dico_ = dico_genename_ens;     
            end
        else  % if model == "Recon3d"
            if data =="Biolinks"
                dico_ = dico_recon3d_ensemble;
            else
                dico_ = dico_recon3d_genename;   
            end
        end    
            
        if HPC == 1
            parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK'))-2) % set the default cores
        end
        
        for c=1:size(discretized,2)
            %clear context_model
            model_files = dir(models_dir);
            file_idx = c;
            file_name= colnames{file_idx};
            file= strcat(string(file_name),'.mat');
            if sum(ismember(string({model_files.name}),file))==0
                % Beging model building for each random class-corhort idx
                [context_model, A] = fastcormics_RNAseq(generic_model, discretized(:,c), rownames,...
                    dico_,biomass_rxn,  already_mapped_tag,consensus_proportion, 1e-4, optional_settings);
                % check model consistency
                %models_keep = zeros(numel(generic_model.rxns), 1); 
                %models_keep(A,1) = 1;
                %context_model = removeRxns(generic_model,generic_model.rxns(setdiff(1:numel(generic_model.rxns),find(models_keep(:,1)))));
                
                % Remove unused genes
                %context_model = removeUnusedGenes(context_model);
                
            % check consistency
            sanity= fastcc_4_fastcormics(context_model,1e-4,0);
            if numel(sanity)==numel(context_model.rxns)
                disp('Consistent Model')
            else
                disp('Inconsistent Model')
            end

            % Adjust Objective function for built models
            if find(ismember(context_model.rxns,biomass_rxn))>0
                model_name = strcat(models_dir,file_name,'.mat');
                context_model = changeObjective(context_model,biomass_rxn,1);
                disp(strcat(model_name,":    ", string(numel(context_model.rxns))))
                if generic_name=="Human1"
                    disp(strcat('Number of sharedEX_to_open rxns is: ',string(numel(intersect(not_medium_constrained, context_model.rxns)))))
                end
            else
                model_name = strcat(models_dir,file_name,'_Inconsistent.mat');
                disp(strcat(model_name,":    ", string(numel(context_model.rxns))))
            end
            
            if HPC==0
                save(model_name,'context_model'); 
            else
                saveModelParfor(model_name,context_model);
            end  
            
            end
        end
        if HPC==1
            delete(gcp); % you have to delete the parallel region after the work is done
        end
    end
end
% 
% 
% model1 = load('models/Human1_TCGA_FPKM_Biolinks_CSF_Thiele2020/TCGA_41_2572_01A_01R_1850_01_Inconsistent.mat');
% model1 = model1.context_model;
% model2 = load('models/Human1_TCGA_FPKM_Biolinks_CSF_Thiele2020/TCGA_41_2572_01A_01R_1850_01_Inconsistent.mat');
% model2 = model2.context_model;
% %not_medium_constrained = {'MAR09153','MAR09167','MAR09269'};%,'MAR09154','MAR10442'};
% intersect(model2.rxns,not_medium_constrained')
% 
% model3 = load('models/Human1_TCGA_FPKM_Biolinks_CSF_Thiele2020/TCGA_41_2572_01A_01R_1850_01_Inconsistent.mat');
% model3 = model3.context_model;
% not_medium_constrained = {'MAR09153';'MAR09167';'MAR09269'};%,'MAR09154','MAR10442'};
% intersect(model3.rxns,not_medium_constrained)
% 
% model4 = load('models/Human1_TCGA_FPKM_Biolinks_CSF_Thiele2020/TCGA_41_2572_01A_01R_1850_01_Inconsistent.mat');
% model4 = model4.context_model;
% not_medium_constrained = {'MAR09153';'MAR09167';'MAR09269';'MAR09154';'MAR10442'};
% intersect(model4.rxns,not_medium_constrained)
% intersect(model4.rxns,'MAR13082')
% intersect(model4.rxns,medium_ex_rxns)
% fastcc_4_fastcormics
% 
% medium_rxns = findRxnsFromMets(model_human1,medium_human1);
% [EX] = findExcRxns(model_human1);
% medium_ex_rxns = intersect(model_human1.rxns(EX),medium_rxns);
% medium_ex_rxns_idx =find(ismember(model_human1.rxns,medium_ex_rxns));
% 
% %medium_not_ex_rxns = setdiff(medium_rxns,model_human1.rxns(EX));
% %medium_not_ex_rxns_idx =find(ismember(model_human1.rxns,medium_not_ex_rxns));
% not_medium_constrained_idx = find(ismember(model_human1.rxns,not_medium_constrained));
% 
% load('matlab.mat', 'model');
% model5 = model;
% idx = find(ismember(model5.rxns,'MAR13082'));
% model5 = changeObjective(model5,model5.rxns(idx));
% solution = optimizeCbModel(model5,'max');
% solution.f
% intersect(model5.rxns,medium_ex_rxns)
% 
% model5.lb(medium_ex_rxns_idx)
% model5.lb(not_medium_constrained_idx)
% EX_names = find(EX);
% y=  model5.lb(setdiff(EX_names,medium_ex_rxns_idx));
% x = EX_names(find(y~=0));
% printRxnFormula(model5,'rxnAbbrList',model5.rxns(x),'metNameFlag',true)
% printRxnFormula(model5,'rxnAbbrList',model5.rxns(find(EX)),'metNameFlag',true)


model4 = load('./Sample_models/Human1_TCGA_FPKM_Biolinks_CSF_Thiele2020_PC/TCGA_06_0747_01A_01R_1849_01.mat');
model4 = model4.context_model;
not_medium_constrained = {'MAR09147';'MAR09153';'MAR09167';'MAR09269';'MAR09276';'MAR09404';'MAR10492';'MAR11895';'MAR13067'};

idx = find(ismember(model4.rxns,'MAR13082'))
idx = find(ismember(model4.rxns,'MAR09167'))
intersect(model4.rxns,not_medium_constrained)