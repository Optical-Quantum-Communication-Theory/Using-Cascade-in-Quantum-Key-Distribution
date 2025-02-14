function extractParameters(paramStruct)
%EXTRACTPARAMETERS Extracts parameters from a parameter struct and places
% them in the workspace where this function is called. Helpful so you can
% use e.g. pz instead of params.pz in many places.
%
% Inputs:
% * paramStruct: the struct of parameters to extract
%
arguments
    paramStruct (1,1) struct
end
    fn = fieldnames(paramStruct);
    for i = 1 : numel(fn)
        fname = string(fn(i));
        field = paramStruct.(fname);
        assignin('caller', fname, field)
    end
end