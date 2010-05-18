function [data_st] = ensemble_filter(data_st,filt)
% Filters data in data_st according to criteria specified in filt.
%
% [data_st] = ensemble_filt(data_st,filt);
%
% Filters data in the data structure (data_st) according to the exclusion and
% inclusion parameters specified in filt.exclude and filt.include.  Each of
% those structures contains either or .and or a .or field or both, which in
% turn contains a set of fields whose names are matched against the
% list of variable names in data_st.vars in order to find the column of data in
% data_st.data to filter.
%
% When specified in the .all structure, all of the conditions in the fields
% must evaluate to true in order for the data to be removed (if part of an
% exclude structure) or retained (if part of an include structure).  When
% specified in a .any structure, any of the conditions have to evaluate to true.
%
% Examples: 
% filt.exclude.any.subject_id = {'tmp_*','01zin79271','08tgs78071'};
% filt.exclude.any.session_id = [1873 1984  1523:1576];
% would cause any rows that have subject IDs beginning with tmp_ to be removed,
% along with the specific subject IDs given by the 2nd and 3rd elements in the
% cell array of strings, as well as any sessions that match the session IDs
% given in the list.
%
% DATE/TIME FILTERING PARAMETERS
% 
% filt.include.all.date_time.start=datenum('01-Jan-1901');
% filt.include.all.date_time.stop=datenum('01-Jan-2020');
% 
% 01/31/07 Petr Janata - adapted from ensemble_apply_crit.m (which didn't have
%                        the added layer of and/or logic)
%
% 02/08/07 Stefan Tomic - added 'exact' argument to strmatch
% 03/15/07 S.T. - added support for using NaN as a filtering criteria
% 09/25/07 PJ - returns if filt spec is empty

%deal with the possibility that params struct was specified
%directly rather than passing "params.filt"
if(isfield(filt,'filt'))
  filt = filt.filt;
end

% If not filtering is specified, then return
if isempty(filt)
  return
end

if isempty(data_st)
  fprintf('%s: empty data struct\n', mfilename);
  return
end

if(iscell(data_st))
  data_st = data_st{1};
end


crit_types_to_proc = fieldnames(filt);
ntypes = length(crit_types_to_proc);

nvars = length(data_st.vars);

% Loop over include and exclude structures
for itype = 1:ntypes
  type_str = crit_types_to_proc{itype};

  % Determine which of the logic operations we're going to perform
  logic_types = fieldnames(filt.(type_str));
  nlog = length(logic_types);
  
  if ~all(ismember(logic_types,{'all','any'}))
    msgstr = sprintf(['ensemble_filter: Found logic types other than ''all'' and' ...
	  ' ''any''\n']);
    error(msgstr)
  end
  
  for ilog = 1:nlog
    logic_str = logic_types{ilog};
    
    % Get a list of the fields to construct masks for
    flds = fieldnames(filt.(type_str).(logic_str));
    nflds = length(flds);

    % Loop over all of the fields associated with this criterion type
    curr_mask = [];
    for ifld = 1:nflds
      fld_str = flds{ifld};

      % Find the field string in the list of variable names
      data_col = strmatch(fld_str,data_st.vars,'exact');
      if isempty(data_col)
	fprintf('ensemble_filter: Did not find criterion field (%s) in list of variables\n',fld_str);
	continue
      end
      
      % Check to see if fld_str is a structure containing limits
      if isstruct(filt.(type_str).(logic_str).(fld_str))
	limit_flds = fieldnames(filt.(type_str).(logic_str).(fld_str));
	nlim = length(limit_flds);
	
	tmp = true(size(data_st.data{data_col}));
	for ilim = 1:nlim
	  limit_str = limit_flds{ilim};
	  crit_val = filt.(type_str).(logic_str).(fld_str).(limit_str);
	  
	  % Make sure there is only one criterion value
	  if length(crit_val) > 1
	    fprintf('ensemble_filter: Too many criterion values\n');
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
	crit_vals = filt.(type_str).(logic_str).(fld_str);
	
	% Make sure some criteria are specified for this field
	if isempty(crit_vals)
	  continue
	end
	
	%if the crit_val is NaN, need to use isnan function,
        %otherwise we can use ismember
	if(~iscell(crit_vals) && any(isnan(crit_vals)))
	  tmp = isnan(data_st.data{data_col});
    else
	  tmp = ismember(data_st.data{data_col}, crit_vals);
	end
	
	% Check to see if any of the criterion values have wildcards, in which case
	% we need to switch to regexp
	if iscellstr(crit_vals) || ischar(crit_vals)
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
    
    %
    % Now construct the final mask that will be used for removing unwanted data
    %
    
    if strcmp(type_str,'include')
      if strcmp(logic_str,'all')
	mask_vect = all(curr_mask,2);
      else
	mask_vect = any(curr_mask,2);
      end
      mask_vect = ~mask_vect; % toggle the mask to be an exclusion mask
    else
      if strcmp(logic_str,'all')
	mask_vect = all(curr_mask,2);
      else
	mask_vect = any(curr_mask,2);
      end
    end
    
    % Perform row extraction  
    for ivar = 1:nvars
      data_st.data{ivar}(mask_vect,:) = [];
    end
  end % for ilog
end % for itype=
