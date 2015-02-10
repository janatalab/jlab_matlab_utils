function sorted_st = ensemble_sort(data_st, sort_var, sort_dir)
% Sorts an Ensemble datastruct 
%
% sorted_st = ensemble_sort(data_st, sort_var, sort_dir)
%
% INPUT:
%   data_st: the Ensemble data structure to be sorted
%   sort_var: the name of the variable to be sorted by
%   sort_dir: the direction ('ascend','descend') to be sorted by. Default
%   is 'ascend'

% 26Jun2014 - Petr Janata
% 09Feb2015 - PJ fixed handling for cells

% Check to make sure that the variable we want to sort by exists
if ~strcmp(sort_var,data_st.vars)
  error('Desired variable for sort (%s) not present in data struct', sort_var)
end

cols = set_var_col_const(data_st.vars);

if ~exist('sort_dir','var') || isempty(sort_dir)
  sort_dir = 'ascend';
end

% Copy the data_st to sorted_st
sorted_st = data_st;
sorted_st.data = cell(size(sorted_st.data)); % Zero out the data portion

% Get the sorted indices
if iscell(data_st.data{cols.(sort_var)})
  [~,sort_idxs] = sort(data_st.data{cols.(sort_var)});
else
  [~,sort_idxs] = sort(data_st.data{cols.(sort_var)},sort_dir);
end

% Loop over all the variables and apply the sort
nvar = length(data_st.vars);
for ivar = 1:nvar
  sorted_st.data{ivar} = data_st.data{ivar}(sort_idxs);
end

return