% parse_brodmann.m
%
% Parses the Brodmann map that is distributed with MRICRO into masks for each
% region. 

% 10/01/07 Petr Janata

srcpath = '/afs/cmb.ucdavis.edu/share/petrlab/mrimasks/';
destpath = srcpath;
srcimg_fname = fullfile(srcpath,'xbrodmann.img');

% Load the image data
V = spm_vol(srcimg_fname);
Y = spm_read_vols(V);

% Get the list of unique values greater than 1
unique_vals = unique(Y(:));

brodmann_areas = unique_vals(unique_vals > 0);

num_areas = length(brodmann_areas);
for iarea = 1:num_areas
  curr_ba = brodmann_areas(iarea);
  
  fprintf('Working on Brodmann Area %d\n', curr_ba);
  
  Yout = zeros(size(Y));
  
  Yout(Y==curr_ba) = 1;
  
  Vout = V;
  Vout.fname = fullfile(destpath,sprintf('BA%d_mask.img', curr_ba));
  Vout.descrip = sprintf('BA%d mask based on Mricro Brodmann map', curr_ba);
  
  fprintf('\tWriting file: %s\n', Vout.fname);
  Vout = spm_write_vol(Vout,Yout);
  
end

