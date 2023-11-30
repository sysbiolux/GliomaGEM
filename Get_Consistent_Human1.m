%Load Human1 model downlaoded on 25/01/2021 from
cd("./GliomGEM")
% https://github.com/SysBioChalmers/Human-GEM/blob/master/model/Human-GEM.mat
load('./Generic_Models/Human-GEM.mat')
dico.genes = ihuman.genes;
ihuman=fixIrr_rFASTCORMICS(ihuman);
[ihuman] = generateRules(ihuman, 0);
% check model consistency
A=fastcc(ihuman,1e-4);
ihuman=removeRxns(ihuman, ihuman.rxns(setdiff(1:numel(ihuman.rxns),A)));
% Remove unused genes
ihuman = removeUnusedGenes(ihuman);
save('./Generic_Models/Human-GEM_Consistent_geneSymbol.mat','ihuman')