function struct1 = struct_union(struct1,struct2)
%
% finds any fields in struct2 that are not in struct1 and copies them to struct1 
%
% The values of the fields that both structs have in common will not be overwritten
% in struct1. The fields missing from struct1 will be populated
% with the values in struct2. This function is useful for
% populating a parameter structure with default values which may be
% missing from the specified parameters.
%
% Copyright (c) 2008 The Regents of the University of California
% All Rights Reserved.
%
% Authors:
% 10/9/2008 - Stefan Tomic, first version.
%


fnamesStruct1 = fieldnames(struct1);
fnamesStruct2 = fieldnames(struct2);

nFields = length(fnamesStruct2);
for iField = 1:nFields
  thisField = fnamesStruct2{iField};
  if(~isfield(struct1,thisField))
    struct1.(thisField) = struct2.(thisField);
  elseif(strcmp(class(struct1.(thisField)),'struct') && ...
	 strcmp(class(struct2.(thisField)),'struct'))
    %if corresponding fields in both structs are themselves structs
    %(substructs), then we'll perform a struct_union on them. If
    %one of them is not a struct (this should be unlikely), then
    %we'll bypass this altogether
    struct1.(thisField) = struct_union(struct1.(thisField),struct2.(thisField));  
  end
end
