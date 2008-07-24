function dater = QA(inDir)
%
% QA(inDir)
%
% Souheil J. Inati  
% Dartmouth College 2000

% Create the name of the first file in inDir
firstfile = fullfile(inDir,'I.001');

% Read the header
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(firstfile);

nX =  im_hdr.imatrix_X; % X Voxels
nY =  im_hdr.imatrix_Y; % Y Voxels
nZ =  im_hdr.slquant;   % Z Voxels
volSize = [nX nY nZ];
adwcount = im_hdr.user9;

% Loop over the volumes until there is no more looping to be done.

% Loop over the passes (volumes)
% Initialize
done = 0;
imageVol = zeros(nX,nY,nZ);
dater = zeros(nX,nY,nZ,adwcount);
for i = 1:adwcount
  % Read in the volume
  [imageVol, lastfile] = GE_readVolume(inDir, i+1, volSize, 16, im_offset);
  dater(:,:,:,i) = imageVol;
  % Done with vol
  fprintf('read %s\n',lastfile);

end % Run loop

return