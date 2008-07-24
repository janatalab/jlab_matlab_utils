function [V] = load_dir(dir_name, filtstr,dims)
% [V] = load_dir(dir_name, filtstr, dims)
%
% Loads .img data from directory in dir_name
%
% Uses strmatch and filtstr to find all .img files that start with filtstr
% dims -- size of image volume (x,y,z).  If this argument is omitted, the
% script defaults to 64x64x27
%

% 05/18/00 PJ

d = dir(dir_name);

%
%  Construct list of .img files
%  Haven't figured out a one-liner to do this
%

fidx = find(~[d.isdir]);

fnames = char(d.name);

img_idx = [];
for f = fidx
  if findstr(fnames(f,:),'.img')
    img_idx(end+1) = f;
  end
end

if nargin > 1
  img_idx = img_idx(strmatch(filtstr,fnames(img_idx,:)));
end

nimg = length(img_idx);
disp(sprintf('Found %d .img files', nimg));

if nargin < 3
  dims = [64 64 27];
end
  
V = zeros(dims(1),dims(2),dims(3),nimg);

disp('Loading files ...')

for iimg = 1:nimg
  curr_fname = fullfile(dir_name, deblank(fnames(img_idx(iimg),:)));
  fid = fopen(curr_fname,'rb','l');
  if fid == -1
    error(sprintf('Could not open file %s', curr_fname))
  end
  V(:,:,:,iimg) = reshape(fread(fid,inf,'integer*2'), size(V,1),size(V,2),size(V,3));
  fclose(fid);
end