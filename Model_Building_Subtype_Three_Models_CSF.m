HPC =0 ; % paramater for running parfor on HPC
datas = ["Rahman2015"];
models = ["Recon3D"];
curation = '_CSF_Thiele2020';

% avoid building model for recon 2 with _DrugMetabolismRemoved since it
% doesn't have this pathway
if contains(curation,"_DrugMetabolismRemoved")==1
    models = ["Recon3D","Human1"];
end

changeCobraSolver('ibm_cplex');
solverOK=changeCobraSolver('ibm_cplex','all');


load('./Generic_models/dico_short.mat')
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
    model_recon3d = load('./Generic_models/Recon3DModel_301_NoDrugMetabolism.mat');
    model_recon3d = model_recon3d.model_recon3d_c;
    model_human1 = load('./Generic_models/Human-GEM_1.5_Consistent_NoDrugMetabolism.mat');
    model_human1 = model_human1.model_human1_c;  
else
    model_recon3d = load('./Generic_models/Recon3DModel_301.mat');
    model_recon3d = model_recon3d.Recon3DModel;
    model_human1 = load('./Generic_models/Human-GEM_1.5_Consistent.mat');
    model_human1 = model_human1.ihuman_consistent;
end

model_recon2 = load('./Generic_models/consistRecon2_4.mat');
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
medium_df = [medium_df;'o2';'co2']
medium_df_recon = strcat(medium_df,'[e]');
medium_recon2 = intersect(medium_df_recon,model_recon2.mets);
medium_recon3d = intersect(medium_df_recon,model_recon3d.mets);

dico_human1_met = readtable('./Generic_models/metabolites_09_2021.tsv','FileType','text');
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
% optional_settings_human1.func = {'c_human','EX_atp[e]'};
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
    optional_settings_human1.medium = medium_human1;
    unpenalizedSystems = {'Transport reactions'};%;'Exchange/demand reactions'};
    x = convertCharsToStrings(model_human1.subSystems);
    x = model_human1.subSystems;
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
    optional_settings_human1.not_medium_constrained = not_medium_constrained;
end
%% Discretization
discretized_biolinks = discretize_FPKM_skewed(fpkm_biolinks,colnames_biolinks);
discretized_rsem = discretize_FPKM(fpkm_rsem,colnames_rsem);
discretized_rahman = discretize_FPKM(fpkm_rahman,colnames_rahman);

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
dico_human = readtable("./Generic_models/genes_08_2021.csv", opts);
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

models_dir = strcat('Sample_models/Subtypes_models/');%,curation,'/');
if ~exist(models_dir, 'dir')
   mkdir(models_dir)
end

meta = readtable('./data/TCGA_TCGBiolinks_metadata_Summary.csv');
meta.sample = replace(meta.sample,'-','_');
meta.sample = extractBefore(meta.sample,16);
Classes = ["AST_IDH_mut",'GBM_IDH_wt','ODG_IDH_mut_Codel','CTRL'];
Classes = ["CTRL"];

for d=1:numel(datas)
    for m=1:numel(models) 
        clear generic_name generic_model data optional_settings biomass_rxn dico_ discretized
        generic_name = models{m};
        data = datas{d};
        
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
           colnames  = extractBefore(colnames_biolinks,16);          
        elseif data== "RSEM"
           discretized = discretized_rsem;
           rownames = rownames_rsem;
           colnames  = colnames_rsem;              
        elseif data== "Rahman2015"
           discretized = discretized_rahman;
           rownames = rownames_rahman;
           colnames  = extractBefore(colnames_rahman,16);  
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
            
%         if HPC == 1
%             parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK'))-2) % set the default cores
%         end
%         
        for i=1:numel(Classes)
            clear context_model
            class = Classes(i);
            %% Model parameters
            if class=="CTRL"
                biomass_rxn = {'biomass_maintenance'};
                optional_settings = optional_settings_recon3d;
                optional_settings = rmfield(optional_settings,'func');
                optional_settings.func = {'biomass_reaction','DM_atp_c_','biomass_maintenance'}; % forced additional reactions into the  model
            else
               optional_settings = optional_settings_recon3d;
               biomass_rxn= biomass_rxn_recon3d;               
            end

            class_idx = find(contains(meta.WHO_2021,class));
            class_barcodes = meta.sample(class_idx)';
            class_idx_2 = find(ismember(colnames,class_barcodes));
            numel(class_idx_2)
            % Begiinng model building for each random class-corhort idx
            [context_model, A] = fastcormics_RNAseq(generic_model, discretized(:,class_idx_2), rownames,...
                dico_,biomass_rxn,  already_mapped_tag,consensus_proportion, 1e-4, optional_settings);
        
            % check consistency
            sanity= fastcc_4_fastcormics(context_model,1e-4,0);
            if numel(sanity)==numel(context_model.rxns)
                disp('Consistent Model')
            else
                disp('Inconsistent Model')
            end
            % Remove unused genes
            context_model = removeUnusedGenes(context_model);

            % Adjust Objective function for built models
            if find(ismember(context_model.rxns,biomass_rxn))>0
                model_name = strcat(models_dir,models{m},"_",datas{d},"_",curation,"_",class,'.mat');
                context_model = changeObjective(context_model,biomass_rxn,1);
                disp(strcat(model_name,":    ", string(numel(context_model.rxns))))
 
            else
                model_name = strcat(models_dir,models{m},"_",datas{d},"_",curation,"_",class,'_Inconsistent.mat');
                disp(strcat(model_name,":    ", string(numel(context_model.rxns))))
            end
            if HPC==0
                save(model_name,'context_model', '-v7' ); 
            else
                %saveModelParfor(model_name,context_model);
            end   
        end
        if HPC==1
            delete(gcp); % you have to delete the parallel region after the work is done
        end
    end
end

%% Mapping back the transcripts to Recon models 
feature astheightlimit 2000
models_dir = 'Sample_models/Subtypes_Models_w_Transcripts/';
if ~exist(models_dir, 'dir')
   mkdir(models_dir)
end

model_recon3d = load('./Generic_models/Recon3DModel_301.mat');
model_recon3d = model_recon3d.Recon3DModel;
model_recon2 = load('./Generic_models/consistRecon2_4.mat');
model_recon2 = model_recon2.model;
model_recon2.description = "recon2_model";

folders = dir('./Sample_models/');
folders = {folders.name};
folders = folders(contains(folders,'Subtypes_Models'));
folders = folders(~contains(folders,'.csv'));

for c=1:numel(folders)
    % Calculate median, minimum and maximum for all models
    model_dir  = string(folders(c));
    model_files = dir('./Sample_models/'+ model_dir+'/');
    model_files = {model_files.name};
    model_files = model_files(contains(model_files,'.mat'));
    model_files = model_files(~contains(model_files,'Inconsistent'));

    for i=1:numel(model_files)
        clear context_model
        model_files{i};
        context_model = load('./Sample_models/'+ model_dir+'/'+string(model_files{i}));
        context_model = context_model.context_model;
        generic_name = split(model_files{i},'_');
        generic_name =generic_name{1};


       if find(context_model.c)>0
            if generic_name=="Human1"
                context_model = changeObjective(context_model,'MAR13082',1);
           else
                context_model = changeObjective(context_model,'biomass_maintenance',1);
            end

            %% Mapping back the transcripts to Recon models 
            if generic_name=="Recon3D"
                context_model_w_transcript = mapModelTranscript(context_model,model_recon3d);

            elseif generic_name=="Recon2"
                context_model_w_transcript = mapModelTranscript(context_model,model_recon2);
            else
              context_model_w_transcript = context_model;
              
            end
       end
       model_file = strcat('./Sample_models/', 'Subtypes_Models_w_Transcripts','/',string(model_files{i}));
       save(model_file,'context_model_w_transcript', '-v7' ); 
    end
end


generic_model = load('./Generic_models/Recon3DModel_301.mat');
generic_model = generic_model.Recon3DModel;
context_model = load('./Sample_models/Subtypes_models/Recon3D_Rahman2015__CSF_Thiele2020_ODG_IDH_mut_Codel.mat');
context_model = context_model.context_model;
context_model_w_transcript = load('./Sample_models/Recon3D_Rahman2015__CSF_Thiele2020_ODG_IDH_mut_Codel.mat');
context_model_w_transcript = context_model_w_transcript.context_model_w_transcript;

table([context_model.genes,context_model_w_transcript.genes])
table([context_model.rxns,context_model_w_transcript.rxns])
numel(context_model.rxns)
numel(context_model_w_transcript.rxns)
