function newstruct = ensemble_tree2datastruct(generic_struct,params)
% Converts a generic data structure to conform to the ensemble data structure
%
% SUPPORTED PARAMS:
%    params.encapsulate_in_cells (0 or 1). Default is 0.
%           This param dictates whether to encapsulate a single vector in a
%           cell. Default is off. This is useful when it is not known
%           ahead of time whether you will be dealing with one or more than
%           one data set in your structure.
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 15, Mar. 2007 - Stefan Tomic
% 5 July, 2007 - Stefan Tomic, edited to conform to current
%                ensemble spec. (data fields do not have to be the same length and
%                their indices do not have to be associated). Consequently,
%                vectors, matrices do not have to be encapsulated in
%                cells. Additionally, an array of generic structs is simply converted to a
%                cell array of ensemble structs
% 26 June, 2008 - Stefan Tomic, arrays of structs now output as a
%                 single Ensemble data struct rather than an array of Ensemble data structs
% 16, Dec. 2009 - Stefan Tomic, bugfix - data are no longer encapsulated in
%                 inner cells, unless you pass an array of structs. e.g.,
%                  if 
%                     generic_struct.fld1 = [1 2 3 4 5]
%                  then
%                    result_struct.vars = {'fld1'}
%                    result_struct.data = {[1 2 3 4 5]}
%		   However if
%		     generic_struct(1).fld1 = [1 2 3 4 5];
%		     generic_struct(2).fld1 = [6 7 8 9 10];
%		   then
%  		     result_struct.vars = {'fld1'}
%		     result_struct.data = { {[1 2 3 4 5],[6 7 8 9 10]} }
%
% Jan 4, 2010 - Stefan Tomic - added encapsulate_in_cells params
%                              Turning this on will encapsulate
%                              any vectors in a cell (even if there is a single
%                              vector). Default is off (do not encapsulate
%                              single vector).
  
  
  
  
try
  params.encapsulate_in_cells;
catch
  params.encapsulate_in_cells = 0;
end
  
%only perform conversion if it looks like this is not an Ensemble
%data struct (i.e. the 'name','type','vars', and 'data' fields are
%not present)
if(all(ismember({'name','type','vars','data'},fieldnames(generic_struct))))
  newstruct = generic_struct;
  return
end

fnames = fieldnames(generic_struct);

typeIdx = strmatch('type',fnames,'exact');
nameIdx = strmatch('name',fnames,'exact');
reportIdx = strmatch('report',fnames,'exact');
metaIdx = strmatch('meta',fnames,'exact');

%if any fieldnames are reserved fieldnames in the datastruct
%(e.g. name or type), assign these to the appropriate places
%instead of to 'vars' and 'data'
if(~isempty(typeIdx))
  fnames = setdiff(fnames,'type');
  type = generic_struct(1).type;
else
  type = '';
end

if(~isempty(nameIdx))
  fnames = setdiff(fnames,'name');
  name = generic_struct(1).name;
else
  name = '';
end

if(~isempty(reportIdx))
  fnames = setdiff(fnames,'report');
  report = generic_struct(1).report;
else
  report = struct();
end

if(~isempty(metaIdx))
  fnames = setdiff(fnames,'meta');
  meta = generic_struct(1).meta;
else
  meta = struct();
end

%if this is an array of structs, we will encapsulate the values corresponding
%to each struct in a cell.
if(length(generic_struct) > 1 || params.encapsulate_in_cells)
  encapsulate = 1;
else
  encapsulate = 0;
end

for ifld = 1:length(fnames)
 
  struct_idx = 1;
  
  for struct_idx = 1:length(generic_struct)
  
    dataCell = generic_struct(struct_idx).(fnames{ifld});

    fieldType = class(dataCell);
    
    if(struct_idx == 1 && ifld == 1)
      newstruct = create_newstruct(fnames,name,type,meta,report);
    end
    
    switch fieldType
	
     case 'struct'
      if(isempty(dataCell))
	assignData = ensemble_init_data_struct;
      else
	assignData = ensemble_tree2datastruct(dataCell);
      end
      
     otherwise
      assignData = dataCell;
      
    end
	
    if(encapsulate)
      newstruct.data{ifld}{struct_idx,1} = assignData;
    else
      newstruct.data{ifld} = assignData;
    end
    
  end

end

return

function newstruct = create_newstruct(fnames,name,type,meta,report)
newstruct = ensemble_init_data_struct;
newstruct.vars = fnames';
newstruct.name = name;
newstruct.type = type;
newstruct.meta = meta;
newstruct.report = report;
return
