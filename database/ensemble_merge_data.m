function outData = ensemble_merge_data(inData,params)
% Merges data from multiple structs into one struct.
%
% outData = ensemble_merge_data(inData,params)
%
%
% accepts a cell array of multiple ensemble data structs, where the
% data structs are all of the same type, with the same number of
% vars and matching var names.
%
% The data from each struct is merged into a single struct and
% returned
%
% Parameters supported:
%    uniqueKeyField: (string) specifies a variable name which
%                             identifies unique row
%                             information. This is useful for merging
%                             multiple data structs with
%                             overlapping information when there
%                             is a key field that describes
%                             repeated row information (e.g. stimulus_id)
%
% 5/24/07 - First Version, S.T.
% 8/13/07 - S.T. added support for uniqueKeyField

% if only a single data struct was passed in, just return the data
% struct with no changes.
if(length(inData) == 1)
  outData = inData;
  return
end



outData = ensemble_init_data_struct;
outData.vars = inData{1}.vars;
outData.data = cell(size(outData.vars));
outVarNames = outData.vars;
outDataCols = set_var_col_const(outData.vars);
outData.type = inData{1}.type;

for iStruct = 1:length(inData)
  
  if(~isempty(setxor(inData{iStruct}.vars,outData.vars)))
    error('All data structs must have matching vars.');
  end
  
  inDataCols = set_var_col_const(inData{iStruct}.vars);
  
  if(isfield(params,'uniqueKeyField'))
    uniqueField = params.uniqueKeyField;
    [keyVals,keyIdxs] = setdiff(inData{iStruct}.data{inDataCols.(uniqueField)},...
				outData.data{outDataCols.(uniqueField)});
  end  
  
  for iData = 1:length(outVarNames)
    
    if(~isfield(params,'uniqueKeyField'))
      nRows = size(inData{iStruct}.data{inDataCols.(outVarNames{iData})},1);
      rowIdxs = 1:nRows;
    else
      nRows = length(keyIdxs);
      rowIdxs = keyIdxs;
    end
    
    %if the data is a row vector, then it will be appended a row at
    %a time, so that we end up with a matrix. If it is a column
    %vector, we'll just append the columns together so they become
    %a long column matrix. This is sort of a cheap and easy way of dealing
    %with whether you want to end up with a matrix or a long
    %vector.
    nCols = size(inData{iStruct}.data{inDataCols.(outVarNames{iData})},2);
    colIdxs = 1:nCols;
    
    if(nRows > 0)
      %vars might not be in the same order in the struct, 
      %so deal with that here
      outData.data{outDataCols.(outVarNames{iData})}(end+1:end+nRows,1:nCols) = ...
	  inData{iStruct}.data{inDataCols.(outVarNames{iData})}(rowIdxs,colIdxs);
    end
    
  end
  
end







