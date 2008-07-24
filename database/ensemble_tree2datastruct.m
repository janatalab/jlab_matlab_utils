function newstruct = ensemble_tree2datastruct(generic_struct,params)
% Converts a generic data structure to conform to the ensemble data structure
%
% params struct isn't used yet and is just there as a placeholder
% for possible future paramters.
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

%if any fieldnames are reserved fieldnames in the datastruct
%(e.g. name or type), assign these to the appropriate places
%instead of to 'vars' and 'data'
if(~isempty(typeIdx))
  fnames = setdiff(fnames,fnames{typeIdx});
  type = generic_struct(1).type;
else
  type = '';
end

if(~isempty(nameIdx))
  fnames = setdiff(fnames,fnames{nameIdx});
  name = generic_struct(1).name;
else
  name = '';
end


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
	if(length(generic_struct) > 1)
	  %multiStructArray is a special case where there is a
          %generic_struct with length >1 and a field which contains
          %another struct array (e.g. an array of form structs,
          %where each form has an array of question structs)
	  fieldType = 'multiStructArray';
	else
	  fieldType = 'struct';
	end
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

    if(struct_idx == 1 && ifld == 1)
      newstruct = create_newstruct(fnames,name,type);
    end
    
    switch fieldType
	
     case 'numeric'
      if(isempty(dataCell))
	assignData = NaN;
      else
	assignData = dataCell;
      end
      newstruct.data{ifld}(struct_idx,1) = assignData;
     
     case 'struct'
      if(isempty(dataCell))
	assignData = ensemble_init_data_struct;
      else
	assignData = ensemble_tree2datastruct(dataCell);
      end
	
      newstruct.data{ifld} = assignData;
      
     case 'multiStructArray'
      if(isempty(dataCell))
	assignData = ensemble_init_data_struct;
      else
	assignData = ensemble_tree2datastruct(dataCell);
      end
      
      newstruct.data{ifld}{struct_idx,1} = assignData;
      
     case 'char'
      if(isempty(dataCell))
	assignData = '';
      else
	assignData = dataCell;
      end
      
      newstruct.data{ifld}{struct_idx,1} = assignData;
      
     case {'cell','matrix'}
      if(isempty(dataCell))
	assignData = [];
      else
	assignData = dataCell;
      end
	
      newstruct.data{ifld}{struct_idx,1} = assignData;
      
    end
	
      
  end
    
  

end

return

function newstruct = create_newstruct(fnames,name,type)
newstruct = ensemble_init_data_struct;
newstruct.vars = fnames';
newstruct.name = name;
newstruct.type = type;
return
