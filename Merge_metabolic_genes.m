changeCobraSolver("ibm_cplex");
model_recon3d = load('./Generic_Models/Recon3DModel_301.mat');
model_recon3d = model_recon3d.Recon3DModel;
model_human1 = load('./Generic_Models/Human-GEM_1.5_Consistent.mat');
model_human1 = model_human1.ihuman_consistent;
model_HMR2 = readCbModel('./Generic_Models/HMR2_MODEL1402200003_url.xml','fileType','SBML');
model_HMR2 = model_HMR2;

load('./Generic_Models/dico_short.mat');
genes_recon3d = model_recon3d.genes;
genes_human1 = model_human1.genes;
genes_hmr2 = model_HMR2.metNames;
genes_hmr2 = genes_hmr2(contains(genes_hmr2,'ENSG'));

genes_recon3d = dico{ismember(dico.ENTREZ,strtok(genes_recon3d,'.')),'SYMBOL'};
genes_human1 = dico{ismember(dico.ENSG,genes_human1),'SYMBOL'};
genes_hmr2 = dico{ismember(dico.ENSG,genes_hmr2),'SYMBOL'};

genes_union = unique([genes_recon3d;genes_human1;genes_hmr2]);

dico_union = dico(ismember(dico.SYMBOL,genes_union),:);
writetable(dico_union,'./Generic_Models/Recon3D_HMR2_Human1_genes.csv','WriteRowNames',true);
