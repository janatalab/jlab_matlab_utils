function [structsAreEqual,firstViolation,violationReason] = compare_structs(struct1,struct2,varargin)
%
% Compares two structs to see if they have identical fields and/or values
%
% struct1 and struct2 are the only required arguments. the rest of the
% arguments are tag,value pairs that specify how the structs will be
% compared.
% 
% supported tags and values:
%    'values' (true or false): specifies whether field values will be
%       compared or not.
%    'substruct' (true or false): specifies whether to check only
%         whether struct1 is contained in struct2 (struct2 simply has
%         more fields and values than struct1 but otherwise are
%         equal). Otherwise the structs will be compared for an exact
%         match.
%    'types' (true or false): specifies whether or not to check the
%         type (class) of the fields for a positive match. This can only
%         be set to false if 'values' is also false. Otherwise, types
%         will automatically be checked.
%    'ignore_fieldnames': a cell array of strings that, if encountered
%         as a fieldname, will be skipped, and will not be compared.
%    'subsetIsOK': a cell array of strings identifying fieldnames for which
%         a subset is of values is OK. Example: If the value of the
%         fieldname 'HalfDecayTimes' in struct1 contains ([0.2 2]), the
%         values of a matched fieldname in struct2 would be deemed equal if
%         it contained the desired values and any others, e.g. [0.2 2 4 8].
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
%                     violation that occurred. This is short and simple to
%                     allow for easy string matching if a more elaborate
%                     report is desired in a higher level script.
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
% Author(s)
% 2007 - Stefan Tomic, first version
% 10/8/2008 - Stefan Tomic, extracted from check_anal_exist and added tag,value functionality
% 12/30/2008 - Fred Barrett, added support for multi-dimensional cell arrays as field values
% 02/24/2009 - Fred Barrett, switched cell checking to compare_cells.m,
%   also added handling of struct arrays
% 11/11/2009 - Fred Barrett, added 'ignore_fieldnames' option
% 01Nov2011 - Petr Janata, added ability to specify those fieldnames for
%             which specification of a subset of values is OK. 
% 30Oct2012 - Improved handling of ignore_fieldnames; eliminated instances
% of strmatch()

numberTypes = {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','double'};
firstViolation = {};
violationReason = '';
VERBOSE = 0;

if (nargin > 2)
  for iarg = 1:2:nargin-2
    switch(varargin{iarg})
     case 'values'
      checkValues = varargin{iarg+1};
     case 'substruct'
      checkSubstruct = varargin{iarg+1};
     case 'types'
      checkTypes = varargin{iarg+1};
     case 'ignore_fieldnames'
			 ignore_fieldnames = varargin{iarg+1};
      case 'subsetIsOK'
        subsetIsOK = varargin{iarg+1};
      case 'verbose'
        VERBOSE = varargin{iarg+1};
    end
  end
end

if ~exist('checkValues','var'), checkValues = true; end
if ~exist('checkSubstruct','var'), checkSubstruct = true; end
if ~exist('checkTypes','var'), checkTypes = true; end
if ~exist('ignore_fieldnames','var')
  ignore_fieldnames = {}; 
elseif ischar(ignore_fieldnames)
  ignore_fieldnames = {ignore_fieldnames};
end
if ~exist('subsetIsOK','var'), subsetIsOK = {}; end

struct1_fieldNames = fieldnames(struct1)';
struct2_fieldNames = fieldnames(struct2)';

struct1_nFields = length(struct1_fieldNames);
struct2_nFields = length(struct2_fieldNames);

fieldNames_intersect = intersect(struct1_fieldNames,struct2_fieldNames);
intersect_nFields = length(fieldNames_intersect);

%the structs must have all fieldnames in common
if(checkSubstruct)
  if(struct1_nFields ~= intersect_nFields)
    fieldNames_diff = setdiff(struct1_fieldNames,struct2_fieldNames);
    for k=1:length(fieldNames_diff)
      if isempty(strmatch(fieldNames_diff{k},ignore_fieldnames))
        firstViolation = fieldNames_diff{k};
        violationReason = 'missing_fields_struct2';
        structsAreEqual = 0;
        return
      end
    end
  end
else
  if(struct1_nFields ~= intersect_nFields || struct2_nFields ~= intersect_nFields)
    fieldNames_xor = setxor(struct1_fieldNames,struct2_fieldNames);
    for k=1:length(fieldNames_xor)
      if isempty(strmatch(fieldNames_xor{k},ignore_fieldnames))
        firstViolation = fieldNames_xor{k};
        violationReason = 'missing_fields_either_struct';
        structsAreEqual = 0;
        return
      end
    end
  end
end

% get the size of the structs ... are they struct arrays?
struct1_size = size(struct1);
struct2_size = size(struct2);

if ~all(struct1_size == struct2_size)
  firstViolation = 'struct2';
  violationReason = 'struct_arrays_different_size';
  structsAreEqual = 0;
  return
else
  nstruct_array = prod(struct1_size); % number of elements
end

% if we are checking only for a substruct, the following will
% still work because we are only cycling through the fields of
% struct1 and checking for a match in struct2. If we are checking 
% for an exact match, then we have already
% verified in the logic above that both structs have the exactly
% the same fields
  
for iField = 1:struct1_nFields
    
  struct1_fieldName = struct1_fieldNames{iField};
  
  if ismember(struct1_fieldName,ignore_fieldnames)
    if VERBOSE
      warning('fieldname %s encountered, but ignored',struct1_fieldName);
    end
    continue
  end
  
  % step through each element of the struct array, compare
  for sIdx = 1:nstruct_array

    struct1_val = struct1(sIdx).(struct1_fieldName);
    struct2_val = struct2(sIdx).(struct1_fieldName);
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
				
				if ~fieldsAreEqual && ismember(struct1_fieldName, subsetIsOK)
					if all(ismember(struct1_val, struct2_val))
						fieldsAreEqual = 1;
					end
				end
 
      case 'char' 
        fieldsAreEqual = strcmp(struct1_val,struct2_val);

      case 'struct'
        [fieldsAreEqual,subFieldViolation,violationReason] = ...
            compare_structs(struct1_val,struct2_val,'values',checkValues,...
            'substruct',checkSubstruct,'types',checkTypes,...
            'ignore_fieldnames',ignore_fieldnames, ...
						'subsetIsOK', subsetIsOK);

      case 'cell'
        [fieldsAreEqual,subFieldViolation,violationReason] = ...
            compare_cells(struct1_val,struct2_val);

      case 'logical'
        if struct1_val == struct2_val
         fieldsAreEqual = 1;
        else
         fieldsAreEqual = 0;
        end

      case 'function_handle'
        struct1_str = func2str(struct1_val);
        struct2_str = func2str(struct1_val);
        fieldsAreEqual = strcmp(struct1_str,struct2_str);
            
      otherwise
        error(sprintf('Fieldtype %s not supported.',struct1_fieldType));


    end % switch struct1

    if (~fieldsAreEqual)
        if(strcmp(struct1_fieldType,'struct'))
          firstViolation = strcat(struct1_fieldName,'.',subFieldViolation);
        else
          firstViolation = {struct1_fieldName};
          violationReason = 'values_differ';
        end
        structsAreEqual = 0;
        return
    end % if (~fieldsAre
  end % for sIdx
  
end % for iField

% if any one field was not equal, then we would not have
% made it this far. Set structsAreEqual =1 and return
structsAreEqual = 1;
return


