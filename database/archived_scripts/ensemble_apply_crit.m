function [data_st] = ensemble_apply_crit(data_st,filt)
% [data_st] = ensemble_apply_crit(data_st,filt);
%
% Filters data in the data structure (data_st) according to the exclusion and
% inclusion parameters specified in filt.exclude and filt.include.  Each of
% those structures contains a set of fields whose names are matched against the
% list of variable names in data_st.vars in order to find the column of data in
% data_st.data to filter.
%
% Example: filt.exclude.subject_id = {'tmp_*','01zin79271','08tgs78071'}
% would cause any rows that have subject IDs beginning with tmp_ to be removed,
% along with the specific subject IDs given by the 2nd and 3rd elements in the
% cell array of strings.
%
% Filtering based on inclusion criteria:  If multiple fields are specified, *all*
% of the conditions have to be true (logical AND) for the data to be retained.
%
% Filtering based on exclusion criteria: Any match on any specified field is
% the basis for exclusion (logical OR).

% 01/25/07 Petr Janata

crit_types_to_proc = fieldnames(filt);
ntypes = length(crit_types_to_proc);

nvars = length(data_st.vars);

for itype = 1:ntypes
  type_str = crit_types_to_proc{itype};
  
  % Get a list of the fields to construct masks for
  flds = fieldnames(filt.(type_str));
  nflds = length(flds);

  % Loop over all of the fields associated with this criterion type
  curr_mask = [];
  for ifld = 1:nflds
    fld_str = flds{ifld};

    % Find the field string in the list of variable names
    data_col = strmatch(fld_str,data_st.vars);
    if isempty(data_col)
      fprintf('Did not find criterion field (%s) in list of variables\n');
      continue
    end
    
    % Check to see if fld_str is a structure containing limits
    if isstruct(filt.(type_str).(fld_str))
      limit_flds = fieldnames(filt.(type_str).(fld_str));
      nlim = length(limit_flds);
      
      tmp = true(size(data_st.data{data_col}));
      for ilim = 1:nlim
	limit_str = limit_flds{ilim};
	crit_val = filt.(type_str).(fld_str).(limit_str);
	
	% Make sure there is only one criterion value
	if length(crit_val) > 1
	  fprintf('ensemble_apply_crit: Too many criterion values\n');
	  continue
	end
	
	% Make sure some criteria are specified for this field
	if isempty(crit_val)
	  continue
	end
	
	switch limit_str
	  case 'start'
	    tmp2 = data_st.data{data_col} > crit_val;
	  case 'start_inc'
	    tmp2 = data_st.data{data_col} >= crit_val;
	  case 'stop'
	    tmp2 = data_st.data{data_col} < crit_val;
	  case 'stop_inc'
	    tmp2 = data_st.data{data_col} <= crit_val;
	end % switch limit_str

	tmp = tmp & tmp2;  % conjoin masks
      end % for ilim
    else
      crit_vals = filt.(type_str).(fld_str);
      
      % Make sure some criteria are specified for this field
      if isempty(crit_vals)
	continue
      end
      
      tmp = ismember(data_st.data{data_col}, crit_vals);
      % Check to see if any of the criterion values have wildcards, in which case
      % we need to switch to regexp
      if iscellstr(crit_vals) | isstr(crit_vals)
	is_wild = ~cellfun('isempty',regexp(crit_vals,'[*]'));
	wild_idxs = find(is_wild);
	for iwild = 1:length(wild_idxs)
	  tmp = tmp|~cellfun('isempty',regexp(data_st.data{data_col}, ...
	      crit_vals{wild_idxs(iwild)})); 
	end
      end
    end % if isstruct(filt.(type_str).(fld_str))
    curr_mask(:,end+1) = tmp;
  end % for ifld
  
  % If it is an inclusion mask, collapse the mask across columns. We need all
  % of the conditions to be true
  if strcmp(type_str,'include')
    mask_vect = all(curr_mask,2);
    mask_vect = ~mask_vect; % toggle the mask to be an exclusion mask
  else
    mask_vect = any(curr_mask,2);
  end
  
  % Perform row extraction  
  for ivar = 1:nvars
    data_st.data{ivar}(mask_vect,:) = [];
  end
  
end % for itype=
