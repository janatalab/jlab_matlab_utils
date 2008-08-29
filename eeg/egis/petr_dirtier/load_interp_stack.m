function [interp_stack,nimg,plot_res]=load_interp_stack(fname, scaling_factor)
% function [interp_stack,nimg,plot_res]=load_interp_stack(fname, scaling_factor);
%
%  interp_stack is a 3D array which is plot_resXplot_resXnum_images
%

% 01/25/00  PJ  modified to return 3D structure
% 06/11/00  PJ  Fixed orientation problem.  Following the last modification,
%               were not being transposed before being placed on the stack.

byte_sex = 'ieee-be';

if nargin < 1
  fid = get_fid('r','*');
  if fid == -1
    disp('Invalid fid in load_interp_stack'),return
  end
  scaling_factor = 500;
else
  fid = fopen(fname,'rb',byte_sex);
  if fid == -1, disp(['Could not open file <' fname '>']),end
end

if nargin < 2, scaling_factor = 500, end

nrow = fread(fid, 1, 'float');
ncolumn = fread(fid,1,'float');
plot_res = fread(fid,1,'float');

% Determine the number of images in the file
fseek(fid,0,'eof');
nbytes = ftell(fid);
nimg = (nbytes-12)/(ncolumn*nrow*4);
fseek(fid,12,'bof');

%
% Original (presumed correct) formulation
%

%for iimg = 1:nimg
%  temp = fread(fid,[ncolumn,nrow],'float');
%  interp_stack(:,(iimg-1)*ncolumn+1:iimg*ncolumn) = temp' / scaling_factor;
%end

interp_stack = zeros(nrow,ncolumn,nimg);

for iimg = 1:nimg
  temp = fread(fid,[ncolumn,nrow],'float');
  interp_stack(:,:,iimg) = temp' / scaling_factor;
end

fclose(fid);

return
