function newstruct = ensemble_tree2datastruct(generic_struct,params)
%
%
% converts a generic data structure to conform to the ensemble data
% structure
%
% params struct isn't used yet. This might be useful for specifying
% the type and names fields of the data struct but 
%
% 15, Mar. 2007 - Stefan Tomic

%only perform conversion if it looks like this is not an Ensemble
%data struct (i.e. the 'name','type','vars', or 'data' fields are
%not present)
if(all(ismember({'name','type','vars','data'},fieldnames(generic_struct))))
  newstruct = generic_struct;
  return
end

fnames = fieldnames(generic_struct);

typeIdx = strmatch('type',fnames,'exact');
nameIdx = strmatch('name',fnames,'exact');
newstruct = ensemble_init_data_struct;

%if any fieldnames are reserved fieldnames in the datastruct
%(e.g. name or type), assign these to the appropriate places
%instead of to 'vars' and 'data'
if(~isempty(typeIdx))
  fnames = setdiff(fnames,fnames{typeIdx});
  newstruct.type = generic_struct(1).type;
end

if(~isempty(nameIdx))
  fnames = setdiff(fnames,fnames{nameIdx});
  newstruct.name = generic_struct(1).name;
end

newstruct.vars = fnames';

for ifld = 1:length(fnames)
 
  fieldType = '';
  struct_idx = 1;
  %have to determine the type of field this is by inspecting the
  %first non-empty value
  while(isempty(fieldType) & struct_idx <= length(generic_struct))
    dataCell = generic_struct(struct_idx).(fnames{ifld});
    if(~isempty(dataCell))
      if(isnumeric(dataCell) & (length(dataCell(:)) == 1))
	fieldType = 'numeric';
      elseif(isnumeric(dataCell) & (length(dataCell(:)) > 1))
	fieldType = 'matrix';
      elseif(isstruct(dataCell))
	fieldType = 'struct';
      elseif(ischar(dataCell))
	fieldType = 'char';
      elseif(iscell(dataCell))
	fieldType = 'cell';
      end
    end
    struct_idx = struct_idx + 1;
  end
  
  %default the fieldType to cell if no data was found
  if(isempty(fieldType))
    fieldType = 'cell';
  end
  
  for struct_idx = 1:length(generic_struct)
  
    dataCell = generic_struct(struct_idx).(fnames{ifld});

    switch fieldType
	
     case 'numeric'
      if(isempty(dataCell))
	assignData = NaN;
      else
	assignData = dataCell;
      end
	
     case 'struct'
      if(isempty(dataCell))
	assignData = {ensemble_init_data_struct};
      else
	assignData = {ensemble_tree2datastruct(dataCell)};
      end
	
     case 'char'
      if(isempty(dataCell))
	assignData = {''};
      else
	assignData = {dataCell};
      end
      
     case {'cell','matrix'}
      if(isempty(dataCell))
	assignData = {[]};
      else
	assignData = {dataCell};
      end
	
    end
	
    newstruct.data{1,ifld}(struct_idx,1) = assignData; 
    
  end
    
  

end

