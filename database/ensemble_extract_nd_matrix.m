function an_st = ensemble_extract_nd_matrix(data_st,params)
% an_st = ensemble_extract_nd_matrix(data_st,params);
%
% Given a list of variables, this function creates an ND matrix, where ND is
% the number of variables with the data in data_st.
%
% A couple of rules are followed. If the data are entirely numeric with a
% single value per variable combination, a numeric ND matrix is returned.
% Otherwise an ND cell array is returned
%
% params should contain the following fields:
% .matrix_dims - specifies which variables in data_st will be used to construct
%                the dimensions of the matrix.
% .dependent_var - specifies which variable provides the values in the matrix.
%                  If this variable is numeric then a numeric matrix is
%                  returned, otherwise a cell array.
%
% Optionally, one can specify:
% .post_process_func - a cell-array of strings that specifies further
% transformations on the matrix before the data are returned.
% .mean_multiple_values - if this is set to 1, and a cell is found to have
% multiple values, the mean of those values will be calculated for that
% cell. The default behavior in this case is to convert the entire data
% matrix into a multi-dimensional cell array.
% 
% NOTE: Currently does not support dynamic handling of response type (radio,
% checkbox, written).  Only handles radio button enums.

% 12/10/08 Petr Janata
% 05/06/09 FB - added params.mean_multiple_values

an_st = ensemble_init_data_struct;

an_st.vars{1} = 'data_matrix';
an_st.vars{2} = 'dim_names';
an_st.vars{3} = 'dim_values';
ancols = set_var_col_const(an_st.vars);

%
% Make sure parameters have been specified
%
if ~all(isfield(params,{'matrix_dims','dependent_var'})) || isempty(params.matrix_dims) ...
      || isempty(params.dependent_var)
  fprintf('%s: No matrix dimension or dependent variables specified\n', mfilename);
  return
end

if iscell(params.dependent_var)
  params.dependent_var = params.dependent_var{1};
end
all_vars = [params.matrix_dims {params.dependent_var}];

%
% See how many variables we want in our output matrix and make sure they all exist
%
have_var_mask = ismember(all_vars, data_st.vars);
if ~all(have_var_mask)
  fprintf('%s: Could not find %d variables in the datastruct: %s\n', ...
      mfilename, sum(~have_var_mask), ...
      cell2str(all_vars(~have_var_mask),','))
  return
end

num_dims = length(params.matrix_dims);
datacols = set_var_col_const(data_st.vars);

%
% Perform any initial filtering
%
if isfield(params,'filt')
  data_st = ensemble_filter(data_st,params.filt);
end

%
% Based on the dependent variable, determine if the output matrix is likely to be a
% numeric or cell array.  
%
make_numeric = 1;

% Are we still numeric
if ~isnumeric(data_st.data{datacols.(params.dependent_var)})
  make_numeric = 0;
end  

%
% Figure out how many unique values there are along each dimension
%
dimension_vals = cell(1,num_dims);
for idim = 1:num_dims
  % Make a copy of the current data
  curr_data = data_st.data{datacols.(params.matrix_dims{idim})};
  
  % Number of unique values along dimension
  dimension_vals{idim} = unique(curr_data);
  
  % Check for and remove NaNs
  if isnumeric(curr_data) && any(isnan(dimension_vals{idim}))
    fprintf('%s: Found NaNs in dimension: \n', mfilename, params.matrix_dims{idim});
    dimension_vals{idim}(isnan(dimension_vals{idim})) = [];
  end
end

%
% Initialize the output matrix
%
if make_numeric
  data_matrix = zeros(cellfun('length', dimension_vals))+NaN;
else
  data_matrix = cell(cellfun('length', dimension_vals));
end

[dim_idxs{1:num_dims}] = deal(0);
curr_dim = 1;
dim_names = params.matrix_dims;

% Figure out what input data we're going to pass in
indata = data_st.data{datacols.(params.dependent_var)};
curr_mask = ones(size(indata));


%
% Now extract the data. This is a bit tricky because the function that does
% this has to be recursive to handle an arbitrary number of dimensions.
% 
data_matrix = burrow(curr_mask, dimension_vals, curr_dim, dim_names, dim_idxs, ...
    indata, data_matrix, data_st, params);

%
% See if we want to do any post-processing
%
try func_list = params.post_process_func; catch func_list = {}; end
if ~isempty(func_list)
  mod_data = data_matrix;
  for ifunc = 1:length(func_list)
    fh = str2func(func_list{ifunc});
    try mod_data = fh(mod_data);
    catch
      fprintf(['%s: Post processing step (%s) failed. Returning data to ' ...
	    'original state\n'], mfilename, func_list{ifunc});
      mod_data = data_matrix;
      break
    end
  end
  data_matrix = mod_data;
end

% Finalize the output structure
an_st.data{ancols.data_matrix} = data_matrix;
an_st.data{ancols.dim_names} = dim_names;
an_st.data{ancols.dim_values} = dimension_vals;

end % function ensemble_extract_nd_matrix

function outdata = burrow(curr_mask, dimension_vals, curr_dim, dim_names, ...
      dim_idxs, indata, outdata, data_st, params)

  datacols = set_var_col_const(data_st.vars);
  curr_data = data_st.data{datacols.(dim_names{curr_dim})};
  try mmv = params.mean_multiple_values; catch mmv = 0; end

  % Get the number of values we have to traverse for the current dimension
  ndim_idxs = length(dimension_vals{curr_dim});

  for iidx = 1:ndim_idxs
    dim_idxs{curr_dim} = iidx;  % to keep track of where we are overall
    
    % Update the mask with the mask for the current dimension's current index
    curr_mask(:,curr_dim) = ismember(curr_data,dimension_vals{curr_dim}(iidx));

    % Check to see if there is a next dim or whether we have descended all the way
    % down
    if curr_dim+1 <= length(dim_idxs)
      outdata = burrow(curr_mask, dimension_vals, curr_dim+1, dim_names, ...
	  dim_idxs, indata, outdata, data_st, params);
    else
      % grab the data
      
      composite_mask = all(curr_mask,2);
      curr_outdata_idx = sub2ind(size(outdata),dim_idxs{:});
      % If we are building a numeric array, make sure we have only a single
      % value, otherwise convert to a cell array
      num_values = sum(composite_mask);
      if isnumeric(outdata) && num_values > 1
        if ~mmv
  	fprintf('%s: Have to convert numeric array to cell array because more than one value was found\n', mfilename)
    outdata = num2cell(outdata);
        end
      end
            % Copy the data
      if num_values > 0
	if ~isnumeric(outdata)
	  outdata{curr_outdata_idx} = indata(composite_mask);
    else
      if length(indata(composite_mask)) > 1 && mmv
    outdata(curr_outdata_idx) = nanmean(indata(composite_mask));
      else
	outdata(curr_outdata_idx) = indata(composite_mask);
      end
	end
      end
    end
  end % for iidx = 1:num_dims

end % burrow
