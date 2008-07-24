function [sub_mask_mtx, subids] = ensemble_make_sub_masks(subid_vect)
% [sub_mask_mtx, subids] = ensemble_make_sub_masks(subid_vect);
%
% Given a vector of subject IDs, a matrix (sub_mask_mtx) is returned in which each
% column contains a logical value indicating the locations within subid_vect that
% correspond to that subject.  The list of subjects is returned in subids, as
% ordered by unique.
%

% 02/01/07 Petr Janata

subids = unique(subid_vect);
nsub = length(subids);
sub_mask_mtx = false(size(subid_vect,1),nsub);

fprintf('Calculating subject masks for %d subjects\n', nsub);
for isub = 1:nsub
  fprintf('.')
  sub_mask_mtx(:,isub) = ismember(subid_vect, subids(isub));
end
fprintf('\n');
