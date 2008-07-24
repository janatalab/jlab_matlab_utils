function data = shuffle_check(epi_dir, file_stub, nslice);
% data = shuffle_check(epi_dir, file_stub, nslice);
%
% Checks for shuffling of slices in EPI images.
%
% You must specify:
%  1) the EPI directory, e.g. '/data2/modfmri/28jan02PJ/epi/run1'
%  2) the root of the name of the EPI files, e.g. 'r28jan02PJ'
%  3) the number of slices, e.g. 27
%
% A color image will be displayed which shows the sum activity in each slice
% through time.  The color of each slice primarily reflects the number of
% voxels in that slice with signal.  Thus, shuffled or duplicated slices appear
% as a vertical line of shuffled colors.
%
% Slice shuffling may result in an erroneous data mask that SPM and most other
% programs depend on to decide which slices to analyze.  Thus, it is best to
% fix these problems (see replace_vols.m)
%
% For scripts that automate this process over multiple subjects and multiple
% runs, see check_slice_intensity.m
%

% v.1.0 06/20/02 Petr Janata
% v.1.1 08/14/02 PJ -- Fixed occassionally incorrect determination of nvol.
%                      Removed requirement to specify nslice


if nargin < 2
  error(['Too few parameters. You must at least specify the EPI directory and root' ...
	' name'])
end

P = spm_get('Files',epi_dir,sprintf('%s*.img', file_stub));
nvol = size(P,1);

fprintf('Analyzing %d volumes beginning with %s in directory %s\n', nvol, file_stub, epi_dir);

if (nargin < 3)
  V = spm_vol(deblank(P(1,:)));
  nslice = V.dim(3);
end

data = zeros(nslice,nvol);
for ivol = 1:nvol
  V = spm_vol(deblank(P(ivol,:)));

  fprintf(' .')
  for islice = 1:nslice
    tmp = spm_slice_vol(V,spm_matrix([0 0 islice]),V.dim(1:2),0);
    data(islice,ivol) = sum(tmp(:));
  end % for islice=
end % for ivol=
    
fprintf('\n')

imagesc(data)
set(gca,'xtick',[0:10:size(data,2)])
xlabel('Volume')
ylabel('Slice')
title(sprintf('%s/%s*.img',epi_dir, file_stub))
