function [structsAreEqual,firstViolation,violationReason] = compare_structs(struct1,struct2,varargin)
%
% Compares two structs to see if they have identical fields and/or values
%
% struct1 and struct2 are the only required arguments
% the rest of the arguments are tag,value pairs that specify how
% the structs will be compared.
% supported tags and values:
%    'values' (true or false): specifies whether field values will be compared or not.
%    'substruct' (true or false): specifies whether to check only
%         whether struct1 is contained in struct2 (struct2 simply has
%         more fields and values than struct1 but otherwise are
%         equal). Otherwise the structs will be compared for an exact
%         match.
%    'types' (true or false): specifies whether or not to check the
%         type (class) of the fields for a positive match. This can only
%         be set to false if 'values' is also false. Otherwise, types
%         will automatically be checked.
%
% defaults if not specified are 'values',true and 'contained',false
%
%  OUTPUTS
%   structsAreEqual - 1 or 0, whether or not struct equality or
%                     substruct is true
%
%   firstViolation  - returns a cell array of fields (or subfields)
%                     that resulted in the first violation of equality (or
%                     substruct-ness). This isn't a complete list
%                     and is only meant to aid in finding
%                     inequalities. 
%
%   violationReason - provides a concise descriptor of the type of
%                     violation that occurred. This is short and simple to allow for
%                     easy string matching if a more elaborate report is desired in a
%                     higher level script.
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
% Author(s)
% 2007 - Stefan Tomic, first version
% 10/8/2008 - Stefan Tomic, extracted from check_anal_exist and added tag,value functionality
%

numberTypes = {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','double'};
firstViolation = {};
violationReason = '';

if (nargin > 2)
  for iarg = 1:nargin-2
    switch(varargin{iarg})
     case 'values'
      checkValues = varargin{iarg+1};
     case 'substruct'
      checkSubstruct = varargin{iarg+1};
     case 'types'
      checkTypes = varargin{iarg+1};
    end
  end
end
 
try
  checkValues;
catch
  checkValues = true;
end
try
  checkSubstruct;
catch
  checkSubstruct=false;
end
try
  checkTypes;
catch
  checkTypes=true;
end

struct1_fieldNames = fieldnames(struct1)';
struct2_fieldNames = fieldnames(struct2)';

struct1_nFields = length(struct1_fieldNames);
struct2_nFields = length(struct2_fieldNames);

fieldNames_intersect = intersect(struct1_fieldNames,struct2_fieldNames);
intersect_nFields = length(fieldNames_intersect);

%the structs must have all fieldnames in common
if(checkSubstruct)
  if(struct1_nFields ~= intersect_nFields)
    firstViolation = setdiff(struct1_fieldNames,struct2_fieldNames);
    violationReason = 'missing_fields_struct2';
    structsAreEqual = 0;
    return
  end
else
  if(struct1_nFields ~= intersect_nFields || struct2_nFields ~= intersect_nFields)
    firstViolation = setxor(struct1_fieldNames,struct2_fieldNames);
    violationReason = 'missing_fields_either_struct';
    structsAreEqual = 0;
    return
  end
end


% if we are checking only for a substruct, the following will
% still work because we are only cycling through the fields of
% struct1 and checking for a match in struct2. If we are checking 
% for an exact match, then we have already
% verified in the logic above that both structs have the exactly
% the same fields
  
for iField = 1:struct1_nFields
    
  struct1_fieldName = struct1_fieldNames{iField};
  struct1_val = struct1.(struct1_fieldName);
  struct2_val = struct2.(struct1_fieldName);
  struct1_fieldType = class(struct1_val);
  struct2_fieldType = class(struct2_val);
  
  %see if the fields are of the same type
  %if both checkTypes and checkValues are turned off, then bypass
  if(~strcmp(struct1_fieldType,struct2_fieldType) &&...
     (checkTypes || checkValues))
    firstViolation = {struct1_fieldName};
    violationReason = 'type_differs';
    structsAreEqual = 0;
    return
  end
  
  %if we are not comparing values and the type of field is anything
  %but a struct, set them to equal and move on to the next
  %field. Otherwise, we will compare their values in the switch statement
  if(~checkValues && ~strcmp(struct1_fieldType,'struct'))
    fieldsAreEqual = 1;
    continue;
  end
  
  %for each case (except for type struct), we need to see if
  %(checkValues) is true. This actually seemed like the most
  %efficient route because a lot of the necessary logic when not comparing
  %values is still contained in the 'for' loop.
  switch struct1_fieldType
      
   case numberTypes
    if(size(struct1_val) == size(struct2_val))
      fieldsAreEqual = all(struct1_val(:) == struct2_val(:));
    else
      fieldsAreEqual = 0; 
    end
 
   case 'char' 
    fieldsAreEqual = strcmp(struct1_val,struct2_val);
  
   case 'struct'
    [fieldsAreEqual,subFieldViolation,violationReason] = compare_structs(struct1_val,struct2_val,'values',checkValues,'substruct',checkSubstruct,'types',checkTypes);
 
   case 'cell'
    %only cell array of strings is supported
    [c,ia,ib] = intersect(struct1_val,struct2_val);
    
    %make sure all elements are in common and in the same order
    if(length(c) == length(struct1_val) ...
       && length(c) == length(struct2_val) ...
       && all(ia == ib))
      fieldsAreEqual = 1;
    else
      fieldsAreEqual = 0;
    end
        
   otherwise
    error(sprintf('Fieldtype %s not supported.',fieldType));
    
  end

  if (~fieldsAreEqual)
    if(strcmp(struct1_fieldType,'struct'))
      firstViolation = strcat(struct1_fieldName,'.',subFieldViolation);
    else
      firstViolation = {struct1_fieldName};
      violationReason = 'values_differ';
    end
    structsAreEqual = 0;
    return
  end
  
  
end

% if any one field was not equal, then we would not have
% made it this far. Set structsAreEqual =1 and return
structsAreEqual = 1;
return

