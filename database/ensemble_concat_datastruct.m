function result_st = ensemble_concat_datastruct(data_st,params)
% Concatenates a cell array of structs into a single struct.
% result_st = ensemble_concat_datastruct(data_st,params);
%
% Concatenates a cell array of data structures in data_st and returns a single
% data structure.  The function checks to make sure that the variables and
% their orders match exactly across the different data structures.
%
% If params.type_as_var is true, then data structures are organized by type,
% such that each type, e.g. response_data becomes a single variable in the
% resulting structure. Note that multiple data structure types are accommodated
% in this way in a single call to ensemble_concat_datastruct.  It is
% recommended that this variable is set to true, but it is currently set to
% false by default for purposes of backward compatibility.

% 05/08/07 Petr Janata
% 03/12/08 Petr Janata

result_st = ensemble_init_data_struct;
try result_st.name=params.outDataName;
catch result_st.name='concat_datastruct';
end

nstruct = length(data_st);

% See if we want to group the input data structures by their type. If we do,
% then the type becomes the type of the returning data structure (result_st).

try type_as_var = params.type_as_var; catch type_as_var = false; end

% Original versions of this script did not operate based on type
if ~type_as_var

  % Copy the variables from the first data struct
  vars = data_st{1}.vars;
  nvars = length(vars);

  result_st.vars = vars;
  result_st.data = data_st{1}.data;

  if nstruct < 2
    return
  end

  for istruct = 2:nstruct
    for ivar = 1:nvars
      if ~strcmp(data_st{istruct}.vars{ivar},vars{ivar})
	fprintf('%s: variable mismatch: Found %s, expected %s\n', mfilename, ...
	    data_st{istruct}.vars{ivar}, vars{ivar});
	result_st = ensemble_init_data_struct;
	return
      else
	result_st.data{ivar} = [result_st.data{ivar}; data_st{istruct}.data{ivar}];
      end
    end % for ivar
  end % for istruct
  
else
  % Get a list of unique datastruct types
  struct_types = cellfun(@getfield,data_st, ...
      repmat({'name'},size(data_st)),'UniformOutput',false);
  
  unique_types = unique(struct_types);
  num_types = length(unique_types);
  
  result_st.vars = unique_types;
  for itype = 1:num_types
    curr_type = unique_types{itype};

    % Initialize a data struct for this particular data type
    curr_st = ensemble_init_data_struct;
    curr_st.name = sprintf('%s_concat', curr_type);
    curr_st.type = curr_type;
    
    % Find all of the instances in the input data structures that match this type
    st_idxs = strmatch(curr_type, struct_types);
    
    % Copy the variable names and data from the first instance
    vars = data_st{st_idxs(1)}.vars;
    nvars = length(vars);
    curr_st.vars = vars;
    curr_st.data = data_st{st_idxs(1)}.data;
    
    % Match the variables and copy the data for the remaining structures
    for istruct = 2:length(st_idxs)
      for ivar = 1:nvars
	if ~strcmp(data_st{st_idxs(istruct)}.vars{ivar}, vars{ivar})
	  fprintf('%s: variable mismatch: Found %s, expected %s\n', mfilename, ...
	      data_st{st_idxs(istruct)}.vars{ivar}, vars{ivar});
	  result_st = ensemble_init_data_struct;
	  return
	else
	  curr_st.data{ivar} = [curr_st.data{ivar}; ...
		data_st{st_idxs(istruct)}.data{ivar}];
	end
      end
    end % for istruct=
    
    result_st.data{itype} = curr_st;
  end % for itype
end % if type_as_var
