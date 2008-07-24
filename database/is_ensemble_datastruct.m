function tf = is_ensemble_datastruct(inStruct)
% Tests a variable to see if it conforms to the ensemble data struct
% specification. Returns 0 if it doesn't conform or 1 if it does.
%
% 4/28/2007 - First Version, Stefan Tomic

if(~isstruct(inStruct))
  tf = 0;
  return
end

emptyStruct = ensemble_init_data_struct;

emptyFieldNames = fieldnames(emptyStruct);
inFieldNames = fieldnames(inStruct);

if(length(emptyFieldNames) ~= length(inFieldNames))
  tf = 0;
  return
end

for ifld = 1:length(emptyFieldNames)

  if (~isfield(inStruct,emptyFieldNames{ifld}))
    tf = 0;
    return
  end
  
  emptyFieldType = class(emptyStruct.(emptyFieldNames{ifld}));
  inFieldType = class(inStruct.(emptyFieldNames{ifld}));
  
  if(~strcmp(emptyFieldType,inFieldType))
    tf = 0;
    return
  end
  
  
end

tf = 1;
return 