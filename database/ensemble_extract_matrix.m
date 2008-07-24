function mtx = ensemble_extract_matrix(data_st,params)
% Extracts the variable specified in params.extract_var from each of the structs in data_st.
% mtx = ensemble_extract_matrix(data_st,params);
%
% Extracts the variable specified in params.extract_var from each of the data
% structures in data_st. If extract_var is not specified, the response_enum
% variable is returned by default.
%
% For example, if the data field in data_st contains a cell array of data
% structures that correspond to question IDs, ensemble_extract_matrix can be
% used to extract the response_enum vector from each of those data structures
% and return it as a single matrix.

% 02/04/07 Petr Janata

mtx = [];

% Verify that we are dealing with a valid data struct
if ~is_ensemble_datastruct(data_st)
  fprintf('%s: Invalid top-level datastruct\n', mfilename);
  return
end

if nargin < 2 || ~isfield(params,'extract_var') || isempty(params.extract_var)
  extract_var = 'response_enum';
else
  extract_var = params.extract_var;
end

% Check to see whether any of the variables in the top-level data struct
% are themselves data structs
num_ds = length(data_st.vars);
data_st_mask = zeros(1,num_ds);
for ids = 1:num_ds
  data_st_mask(ids) = is_ensemble_datastruct(data_st.data{ids});
end

if ~any(data_st_mask)
  ds_cols = set_var_col_const(data_st.vars);
  mtx = data_st.data{ds_cols.(extract_var)};
  return
else
  for ids = 1:num_ds
    ds_cols = set_var_col_const(data_st.data{ids}.vars);

    tmp = data_st.data{ids}.data{ds_cols.(extract_var)};

    % Initialize the output variable
    if ids == 1
      nrows = length(tmp);
      if iscell(tmp)
        mtx = cell(nrows,num_ds);
      else
        mtx = zeros(nrows,num_ds);
      end
    end

    % Check to make sure that the length of the temporary vector matches the
    % number of rows in the output matrix
    if length(tmp) ~= nrows
      fprintf('%s: number of rows mismatch: Expected %d, found %d\n', mfilename, nrows, length(tmp));
      mtx = [];
      return
    end

    % Copy the source data into the output matrix
    mtx(:,ids) = tmp;

  end % for ids
end % if ~any(data_st_mask)
return

