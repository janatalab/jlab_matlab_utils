function [mask_mtx, unique_vals] = make_mask_mtx(data_vect)
% [mask_mtx, unique_vals] = make_mask_mtx(data_vect);
%
% Given a data vector (data_vect) containing some set of values or identifiers, a
% matrix (mask_mtx) and list of unique values (unique_vals) are returned in which
% each column of mask_mtx contains a logical value indicating the locations within data_vect
% that correspond to a unique values in unique_vals.  The list of unique
% valuess is returned as dered by unique.
%

% 02/14/07 Petr Janata - generalized from ensemble_make_sub_masks

unique_vals = unique(data_vect);
nvals = length(unique_vals);
mask_mtx = false(size(data_vect,1),nvals);

fprintf('Calculating masks for %d unique values\n', nvals);
for ival = 1:nvals
  fprintf('.')
  mask_mtx(:,ival) = ismember(data_vect, unique_vals(ival));
end
fprintf('\n');
