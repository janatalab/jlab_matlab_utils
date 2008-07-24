function myNewStruct = flatten_children_structs(myStruct)
%
% myNewStruct = flatten_children_structs(myStruct)
% 
% accepts a struct with one or more struct arrays in some its
% fields (child struct arrays) and brings the fields of the "child
% structs" up to the top level. A struct array with the same length
% as the child structs is returned. All fields with struct arrays
% must have the same number of elements.
% 
% e.g. if myStruct has the following fields:
%
%      field1: 34
%      field2: [2 3 4 5 6]
%      field3: {'one' 'two' 'three'}
%      field4: [1x300 struct] with fields: field5, field6, field7
%
% the returned myNewStruct is a 1x300 struct with fields:
%
%      field1: 34 (across all elements of array)
%      field2  [2 3 4 5 6] (across all elements of array)
%      field3  {'one' 'two' 'three'} (across all elements of array)
%      field5: <value from field4 struct>
%      field6: <value from field4 struct>
%      field7: <value from field4 struct>
%
% Note that if there are child structs that are not in an array
% these are simply treated as a struct array of length 1 and these
% are flattened into the parent struct.
%
% All child struct arrays must be of the same length.
%
% 26 Jan, 2007 - First version, S.T. 


if(~isstruct(myStruct) | length(myStruct) ~= 1 | ndims(myStruct) > 2)
  error('input argument must be a struct of length 1');
end

%obtain the fieldnames and values of myStruct
myStructFieldNames = fieldnames(myStruct);
myStructFieldValues = struct2cell(myStruct);

%go through all fields and find out which ones are structs
%for those that are structs, record the lengths of the struct arrays
for iField = 1:length(myStructFieldNames)
  fieldIsStruct(iField) = isstruct(myStruct.(myStructFieldNames{iField}));
  if(fieldIsStruct(iField))
    fieldStructLength(iField) = length(myStruct.(myStructFieldNames{iField}));
  else
    fieldStructLength(iField) = 0;
  end
end

%the lengths of all child structs must be equal
if(any(diff(fieldStructLength(fieldStructLength ~= 0)) ~= 0))
  error('Child struct arrays must all have the same length');
else
  childrenLengths = fieldStructLength(find(fieldStructLength ~= 0,1));
end

%obtain the indexes of the fields for the children struct arrays
childStructIdxs = find(fieldIsStruct);

%obtain the fieldnames that are not structs
myNewStructFieldNames = myStructFieldNames(~fieldIsStruct);


for idxChildStruct = 1:childrenLengths
  %obtain the values for fields that are not structs
  myNewStructFieldValues = myStructFieldValues(~fieldIsStruct);

  %construct the new struct array, going through each element of
  %the child struct arrays
  for ichildStruct = childStructIdxs
    childStructName = myStructFieldNames{ichildStruct};
    tmpChildStruct = myStruct.(childStructName);
    childStructFieldNames = fieldnames(tmpChildStruct);
  
    childStructFieldValues = struct2cell(tmpChildStruct(idxChildStruct));
    
    
    
    
    
    if(idxChildStruct == 1)
      myNewStructFieldNames  = {myNewStructFieldNames{:} childStructFieldNames{:}};
    end
    
    
    myNewStructFieldValues = {myNewStructFieldValues{:} childStructFieldValues{:}};
    
    
  end
  
  newStructParams = cell(1,length(myNewStructFieldNames)*2);
  newStructParams(1:2:end) = myNewStructFieldNames;
  newStructParams(2:2:end) = myNewStructFieldValues;
  myNewStruct(idxChildStruct) = cell2struct(myNewStructFieldValues,myNewStructFieldNames,2);
end

return
