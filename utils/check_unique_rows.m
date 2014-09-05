function [non_unique_mask, non_unique_vals] = check_unique_rows(mtx,verbose)
% Checks whether all rows in the input matrix, mtx, are unique, and returns
% a mask of non-unique values and an extracted set of the non-unique values
%
% This function currently operates only on categorical class data
%
% [non_unique_mask, non_unique_vals] = check_unique_rows(mtx);
%
% TODO: Add support for non-categorical variables

% 04Sep2014 - Petr Janata

if nargin < 2
  verbose = 1;
end

if ~any(strcmp({'nominal','ordinal'}, class(mtx)))
  error('Class of input data must be nominal or ordinal')
end

nrows = size(mtx,1);

non_unique_mask = false(nrows,1);
non_unique_vals = [];

uniqueRows = unique(mtx,'rows');
numNonUnique = nrows - size(uniqueRows,1);
if numNonUnique
  if verbose
    fprintf('Found %d non-unique rows%s\n',numNonUnique);
  end
  % Find the non-unique rows
  [~,idxs] = ismember(mtx,uniqueRows,'rows');
  
  % Tabulate idxs to see which appears more than once
  t = tabulate(idxs);
  
  % Display the nonUnique rows
  non_unique_vals = uniqueRows(t(:,2)>1,:);
  if verbose
    non_unique_vals
  end
  non_unique_mask = ismember(mtx,non_unique_vals,'rows');
end

end % function