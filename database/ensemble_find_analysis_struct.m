function an_st_idxs = ensemble_find_analysis_struct(an_st_array, search_crit)
% Returns indices of structs within a cell array that match certain criteria.
%
% an_st_idx = ensemble_find_analysis_struct(an_st_array, search_crit);
%
% Retrieves the indices of analysis structures within a cell array of analysis
% structures that match certain criteria, such as analyis type or values in
% fields of structures contained in the meta field associated with the
% analysis.
%
% This function is very useful for post-processing analyses that have been returned
% by some of the generic lower-level analysis functions such as
% ensemble_enum_hist.
%
% search_crit is a structure in which the field names should refer to fields in
% the data structure (an_st_array), e.g. type, vars, data.  The value(s) of the
% criterion field in the search crit struct are the values of the criterion
% variable to match.  
%
% Example, if you are looking for an analysis structure that has subject_id as
% one of the variables, the criterion field would be 'vars' and its value would
% be 'subject_id'

% 01/31/07 Petr Janata
an_st_idxs = []; % Initialize the output variable to empty

if isstruct(an_st_array)
  an_st_array = {an_st_array};
end

na = length(an_st_array); % number of analyses in the list

% Set up a global mask for keeping track of which analysis structures survive
% all of the criteria
global_mask = true(1,na);

% Get a list of criterion fields that we are going to search on
crit_fld_list = fieldnames(search_crit);
ncrit_fld = length(crit_fld_list);

% Loop over top-level criterion fields
for ifld = 1:ncrit_fld
  crit_fld = crit_fld_list{ifld};
  
  % Get a list of the analysis structures that contain this criterion field
  crit_fld_mask = cellfun(@isfield,an_st_array, repmat({crit_fld},1,na));
  
  % Make sure that at least some analysis structure had this mask
  if ~any(crit_fld_mask)
    fprintf('None of the analysis structures had the criterion field: %s\n', crit_fld);
    an_st_idxs = [];
    return
  end
  
  % Check to see if the criterion field is itself a structure or if it contains
  % a criterion value
  crit_val = search_crit.(crit_fld);
  crit_val_is_struct = isstruct(crit_val);

  % Either way, we want to extract the corresponding values from the analysis structure
  if isstr(crit_val) || crit_val_is_struct
    analysis_vals = cell(1,na);
  else
    analysis_vals = zeros(1,na)*NaN;  % 05/01/07 - this might be a bug to size
                                      % it to na ?
  end
  if isstr(crit_val) || crit_val_is_struct
    analysis_vals(crit_fld_mask) = ...
	cellfun(@getfield,an_st_array(crit_fld_mask), ...
	repmat({crit_fld},1,sum(crit_fld_mask)), 'UniformOutput',false);
  else
    analysis_vals(crit_fld_mask) = ...
	cellfun(@getfield,an_st_array(crit_fld_mask), ...
	repmat({crit_fld},1,sum(crit_fld_maks)));
  end
  
  if isstr(crit_val)
    analysis_vals(~crit_fld_mask) = repmat({''},1,sum(~crit_fld_mask));
  end
  
  % If the criterion value is not a structure, then match the analysis values
  % against the criterion values to get our mask
  curr_mask = false(1,na);
  if ~crit_val_is_struct
    if isstr(crit_val)
      curr_mask = ismember(analysis_vals, crit_val);
    else
      curr_mask = analysis_vals == crit_val;
    end
  else
    % Otherwise, recursively call this function
    good_idxs = ensemble_find_analysis_struct(analysis_vals,crit_val);
    curr_mask(good_idxs) = true;
  end
  
  global_mask = global_mask & curr_mask;
  
  % Terminate the search if there are no more analyses to check out
  if ~any(global_mask)
    an_st_idxs = [];
    return
  end
end % for ifld = 1:ncrit_fld

an_st_idxs = find(global_mask);

return
