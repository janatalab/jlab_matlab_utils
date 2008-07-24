function simple_replace(epi_dir,epi_stub,bad_vol,good_vols)
% simple_replace(epi_dir,epi_stub,bad_vol,good_vols);
%
% Replaces specified volume with average of other specified volumes.
% Good for fixing slice-shuffled volumes.
%
% Assumes that the last portion of the file names has the format:
%     i%04d.img, e.g. i0037.img
%
% WARNING: Make sure you have a copy of the original data, as this will
% overwrite the existing image.

% v.0.0.1 6/20/02 Petr Janata
% v.0.0.2 08/09/02 PJ  Data are read in via SPM functions to avoid potential scaling
%                      factor problems.

if nargin < 4
  error('Too few arguments')
end

nrep = length(good_vols);

bad_vol_fname = spm_get('Files',epi_dir,sprintf('%s*i%04d.img', epi_stub,bad_vol));
if isempty(bad_vol_fname)
  error(sprintf('Could not locate a file in %s starting with %s and volume number %d', epi_dir, epi_stub, bad_vol))
end

Vbad = spm_vol(bad_vol_fname);

good_data = zeros(Vbad.dim(1:3));
nslice = Vbad.dim(3);

for irep = 1:nrep
  good_vol_fname = spm_get('Files',epi_dir,sprintf('%s*i%04d.img',epi_stub,good_vols(irep)));

  if isempty(good_vol_fname)
    error(sprintf('Could not locate a file in %s starting with %s and volume number %d', epi_dir, epi_stub, good_vols(irep)))
  end
  
  Vgood = spm_vol(good_vol_fname);
  for islice = 1:nslice
    good_data(:,:,islice) = good_data(:,:,islice) + ...
	spm_slice_vol(Vgood,spm_matrix([0 0 islice]),Vgood.dim(1:2),0);
  end
end

good_data = good_data/nrep; % get the mean

% Write out the modified volume
Vbad = spm_write_vol(Vbad,good_data);

