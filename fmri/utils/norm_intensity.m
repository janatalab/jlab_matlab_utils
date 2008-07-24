function [nV] = norm_intensity(V)
% [nV] = norm_intensity(V);
%
% Normalizes intensity of each image volume in a series of image volumes
%

% 05/18/00 PJ

thresh_pct = 0.4;

%
%  Set a threshold for the brain mask to use in the intensity normalization
%

dims = size(V);
rV = reshape(V,prod(dims(1:3)),dims(4));

% Get the mean intensity of each image volume
mi = mean(rV);

% Set the threshold for each image volume
thresh = mi*thresh_pct;

% Find all voxels exceeding the threshold
mask = rV >= repmat(thresh,size(rV,1),1);

% Determine the scaling factor for each volume
S = sum(mask)*100./sum(rV .* mask);

% Multiply the entire by the scaling factor
mask_idx = find(mask);

rV = rV .* repmat(S,size(rV,1),1);

% Reshape back to the original size
nV = reshape(rV, dims);

