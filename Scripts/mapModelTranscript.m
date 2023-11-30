function [context_model_w_transcript] = mapModelTranscript(context_model,generic_model)
%mapModelTranscript mapps back transcripts from .genes using a generic
%model
    rxns = context_model.rxns;
    context_model_w_transcript = context_model;
    for i=1:numel(rxns)
        rxns_b_idxs = find(ismember(generic_model.rxns,rxns(i)));
        [~,genes_b_idxs] = find(generic_model.rxnGeneMat(rxns_b_idxs,:));
        [~,genes_a_idxs] = find(context_model.rxnGeneMat(i,:));
        if numel(genes_b_idxs) >0
            genes_b = cellstr(generic_model.genes(genes_b_idxs));
            context_model_w_transcript.genes(genes_a_idxs) = genes_b;

        end
    end
    context_model_w_transcript.genes = cellstr(context_model_w_transcript.genes);
end