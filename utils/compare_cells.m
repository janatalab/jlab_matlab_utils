function [cellsAreEqual,firstViolation,violationReason] = compare_cells(cell1,cell2)

% this function compares the structure and contents of multi-dimensional and nested cells
% 
% this function has been modeled after compare_structs.m (S. Tomic)
% 
% cell1 and cell2 are the only required arguments
% 
% the rest of the arguments are maintained for compatibility with
% compar_struct, and are passed onto compare_struct in the event that cell
% contents are of class 'struct'. tag,value pairs that specify how the
% structs will be compared.
% supported tags and values:
% 
%    'values' (true or false): specifies whether field values will be
%       compared or not.
%    'substruct' (true or false): specifies whether to check only whether
%       struct1 is contained in struct2 (struct2 simply has more fields and
%       values than struct1 but otherwise are equal). Otherwise the structs
%       will be compared for an exact match.
%    'types' (true or false): specifies whether or not to check the type
%       (class) of the fields for a positive match. This can only be set to
%       false if 'values' is also false. Otherwise, types will
%       automatically be checked.
% 
% this function first checks to see if both inputs are cells, and then
% compares the structure (# dimensions and size in each dimension) of the
% provided cells. Any inequality in # of dimensions or size in any
% dimension will be viewed as a cell inequality. If cells are of identical
% structure, they will be converted to a vector [new_vect = cellN(:)], and
% then each element from both cells will be compared on variable class and
% then content.
% 
%  OUTPUTS
%   cellsAreEqual - 1 or 0, whether or not cell equality is true
%   firstViolation - the index in the vector-transformed cells where the
%       first inequality of content or variable class was found. In the
%       case of cells within cells, this function is called recursively.
%       when recursing on self, upon finding an error, the index in the
%       topmost parent cell where the inequality was found is returned in
%       firstViolation. In the case that the structure of the cells
%       differs, or the case that either or both of the inputs are not
%       cells, firstViolation = 0;
%   violationReason - provides a concise descriptor of the type of
%       violation that occurred.
% 
% NOTE: if one cell is a transposition of the other, they will be found
% unequal, either by size, or by contents
% 
% Copyright (c) 2009-2012 The Regents of the University of California
% All Rights Reserved.
%
% FB 2009.02.24

% 
% initialize variables
% 
cellsAreEqual = 0;
firstViolation = 0;
violationReason = '';

numberTypes = {'int8','uint8','int16','uint16','int32','uint32','int64',...
    'uint64','double'};

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

if ~exist('checkValues','var'),    checkValues = true; end
if ~exist('checkSubstruct','var'), checkSubstruct = true; end
if ~exist('checkTypes','var'),     checkTypes = true; end

% 
% check to see if inputs are cells
% 
if ~iscell(cell1) && ~iscell(cell2)
  violationReason = 'neither input is of type ''cell''';
  return
elseif ~iscell(cell1)
  violationReason = 'cell1 is not of type ''cell''';
  return
elseif ~iscell(cell2)
  violationReason = 'cell2 is not of type ''cell''';
  return
end

% 
% compare sizes of cells
% 
s1 = size(cell1);
s2 = size(cell2);

% if cells are not the same exact size, they will be treated as unequal
% NOTE: transposed cell arrays will not be treated as equal
if length(s1) ~= length(s2)
  violationReason = 'cells have different number of dimensions';
  return
end

sm = (s1 == s2);
if ~all(sm)
  violationReason = 'cells have different sizes in one or more dimension';
  return
end

% convert cells to vectors
v1 = cell1(:);
v2 = cell2(:);
nv1 = length(v1);
nv2 = length(v2);
if nv1 ~= nv1
  violationReason = 'vectors (converted cells) have different sizes';
  return
end

% step through cells, compare at each index
for i=1:nv1
    
  % compare cell class
  c1 = class(v1{i});
  c2 = class(v2{i});  
  if isempty(strmatch(c1,c2,'exact'))
    firstViolation = i;
    violationReason = sprintf(['cell contents have different class at '...
        'index %d'],i);
    return
  end

  % compare cell contents
  same = 1;
  switch c1     
      case numberTypes
        if any(size(v1{i}) ~= size(v2{i})) || ~all(v1{i} == v2{i})
          same = 0;
          subVioReason = sprintf('%s differed',c1);
        end
      case 'char'
        if ~strcmp(v1{i},v2{i})
          same = 0;
          subVioReason = 'strings differed';
        end
      case 'struct'
        [structsAreEqual,subStructVio,subVioReason] = compare_structs(...
            v1{i},v2{i},'values',checkValues,'substruct',checkSubstruct,...
            'types',checkTypes);
        if ~structsAreEqual, same = 0; end
      case 'cell'
        [subCellsEqual,subCellVio,subVioReason] = compare_cells(v1{i},v2{i});
        if ~subCellsEqual, same = 0; end
      case 'logical'
        if v1{i} ~= v2{i}
          same = 0;
          subVioReason = sprintf('logical vars differed');
        end          
      otherwise
        same = 0;
        subVioReason = sprintf('class %s not supported',c1);
  end

  if ~same
    firstViolation = i;
    violationReason = sprintf('cell contents differed (%s)',subVioReason);
    return
  end
end % for i=1:nv1

% no differences found, return equality
cellsAreEqual = 1;
