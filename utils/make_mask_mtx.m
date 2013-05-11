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
% 10May2013 PJ - optimized by using == and strcmp instead of ismember

if size(data_vect,2) > 1
	nrows = size(data_vect,1);
	if iscell(data_vect)
		tmpdata = data_vect;
		nvars = size(tmpdata,2);
		numCombos = 0;
		while size(tmpdata,1)
			numCombos = numCombos+1;
			if mod(numCombos,20) == 0
				fprintf('%d', numCombos);
			else
				fprintf('.');
			end
			
			% Get a combination
			currCombo = tmpdata(1,:);
			
			% Store it for future use
			[unique_vals(numCombos,1:nvars)] = deal(currCombo);
			
			mask_mtx(:,numCombos) = all(strcmp(data_vect,repmat(currCombo,nrows,1)),2);
						
			currComboMask = all(strcmp(tmpdata,repmat(currCombo,size(tmpdata,1),1)),2);
			% Eliminate those rows
			tmpdata(currComboMask,:) = [];
		end
		
		fprintf('\n');
		return
	else
		unique_vals = unique(data_vect,'rows');
	end
else
	unique_vals = unique(data_vect);
end

% Eliminate NaNs if the exist
if isnumeric(unique_vals)
  unique_vals(isnan(unique_vals)) = [];
end

nvals = length(unique_vals);
mask_mtx = false(size(data_vect,1),nvals);

fprintf('Calculating masks for %d unique values\n', nvals);
for ival = 1:nvals
	if mod(ival,20) == 0
		fprintf('%d', ival);
	else
		fprintf('.');
  end
  if isnumeric(unique_vals)
    mask_mtx(:,ival) = data_vect == unique_vals(ival);
  elseif iscell(unique_vals(ival)) && ischar(unique_vals{ival})
    mask_mtx(:,ival) = strcmp(unique_vals{ival}, data_vect);
  else
    mask_mtx(:,ival) = ismember(data_vect, unique_vals(ival));
  end
end
fprintf('\n');

return

