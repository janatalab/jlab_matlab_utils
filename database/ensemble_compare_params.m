function found = ensemble_compare_params(filename,params)
% Compares a saved params struct with a given params struct in memory.
%
% Compares a set of params with a params struct saved in the 'data'
% field of a data structure stored in a file. It is assumed that
% the name of the data structure is the same as the params
% subfield.
%
% For example, if a set of ani parameters are passed in as
% params.ani (where params.ani.downSampleFactor and
% params.ani.FirstCBU specify two parameters of the ani), this
% function will attempt to load an 'ani' struct from the file
% specified by the input 'filename' and the 'params' data will be
% compared with the input params.
%
% If params contains any other fields, such as params.ci,
% params.pp, the file will also be searched for ci and pp data
% structures, so it is best to only pass the param fields that you
% wish to be compared with this file.

if(~exist(filename,'file'))
  error(sprintf('Filename %s does not exist.',filename));
end

dataStructNames = fieldnames(params);

for iDataStruct = 1:length(dataStructNames)

  thisDataStructName = dataStructNames{iDataStruct};
  
  fileData = load(filename,thisDataStructName);
  
  thisDataStruct = fileData.(thisDataStructName);
  thisDataStructCols = set_var_col_const(thisDataStruct.vars);
  
  if(~ismember('params',thisDataStructCols))
    disp([sprintf('Filename %s, data struct %s does not contain a',filename,thisDataStructName) ...
		  ' params struct']);
    found = 0;
    return
  end
  
  loadedParams = thisDataStruct.data{thisDataStructCols.params};
  paramFields = fieldnames(params.(thisDataStructName));
  
  
  
end
