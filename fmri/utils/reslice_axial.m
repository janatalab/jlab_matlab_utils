function reslice_axial(fname,VOX,BB)
%
% function reslice_axial(fname,VOX,BB)
%
% Takes a file name fname and reslices it axially using trilinear
% interpolation to the space defined by it. 
% If VOX is not specified it defaults to [1 1 1].
%
% VOX = [3 3 3] for example asks for 3x3x3 mm voxels
% BB = bounding box of new image.  If this is left empty it will use the
% bounding box of the original image.
%
% Uses the spm function spm_write_sn and trilinear interpolation
%
% Souheil Inati
% Dartmouth College

% Modification history:
% 10/12/01 PJ -- Added BB parameter

if (nargin<2)
  VOX = [1 1 1];
end

if nargin<3
  BB = [];
end

V = spm_vol(fname);

% Compute the bounding box
% Find the mins and maxes
if isempty(BB)
  lims1 = V.mat*[1 1 1 1]';
  lims2 = V.mat*[V.dim(1:3) 1]';
  coordlims = [ lims1(1:3) lims2(1:3) ];
  mins = [min(coordlims,[],2)];
  maxs = [max(coordlims,[],2)];
  BB = [mins';maxs'];
end

% Define the identity sn3d matrix and save it next to the image
[imgpath,imgname] = fileparts(fname);
sn3dfilename = fullfile(imgpath,'identity_sn3d.mat');
mgc    = 960209;
MF     = eye(4);
MG     = eye(4);
Affine = eye(4);
Dims   = [1 1 1; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1];
Transform = [];
save(sn3dfilename,'mgc','Dims','Affine','MF','MG','Transform');

% Use Trilinear interpolation
Hold = 1;

% Write it out
spm_write_sn(fname,sn3dfilename,BB,VOX,Hold);

% Clean Up
delete(sn3dfilename)

return
