function match = ensemble_verify_row_match(data_st,params)
% Verifies that rows match across subordinate structs within a given struct.
% match = ensemble_verify_row_match(data_st,params)
%
% Makes sure that rows match across subordinate data structures in the
% top-level data structure (data_st) for all of the variables specified in the
% match_vars fields of the params structure.  If no match_vars field is
% specified, all common variables are matched.
%
% This function is a mechanism to make sure that data are not mis-aligned
% across things like question IDs when the intention is to extract a matrix of
% data from individual children data structures within a parent data structure.

% 02/04/07 Petr Janata

match = true;  % start out optimistic

num_ds = length(data_st.vars);

if nargin < 2 || ~isfield(params,'match_vars') || isempty(params.match_vars)
  match_vars = data_st.data{1}.vars;
else
  match_vars = params.match_vars;
end
nmatch = length(match_vars);

% Find the common set of variables to match on

% Compare the data
for ids = 1:num_ds
  child_cols = set_var_col_const(data_st.data{ids}.vars);
  
  for im = 1:nmatch
    curr_var = match_vars{im};
    
    ep.extract_var = curr_var;
    mtx = ensemble_extract_matrix(data_st,ep);
    
    if isempty(mtx)
      match = false;
      return
    end
    
    if iscellstr(mtx)
      ncol = size(mtx,2);
      for icol = 2:ncol
	if ~all(strcmp(mtx(:,icol),mtx(:,1)))
	  match = false;
	  return
	end
      end
    else
      if any(any(diff(mtx,1,2),2))
	match = false;
	return
      end
    end
  end  
end % for ids

return

