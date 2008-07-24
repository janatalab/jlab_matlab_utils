function [outVol] = GE_convertT2PD(inDir,outstub)
%
% [vols] = GE_convertT2PD(inDir,outstub)
%
% Converts a series of GE slices obtained using a dual-echo protocol into two
% Analyze format files. The first file contains the T2-weighted image, and the
% other contains the proton density (PD) image.  The images are interleaved
% during acquisition.
%
% 09/03/03 Petr Janata
%

% This script was adapted from GE_convertVolume.m:
%
% Converts a series of GE slices into Analyze format for SPM.
% inDir is the name of the directory containing the first file of
% the series, e.g. 003
% runnum is the run number of the set you want to convert
% outName is used for naming the Analyze files
% status = 1 if there is an error.
%
% Souheil J. Inati  
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

img_types = {'PD','T2'};
ntypes = length(img_types);

status = 1;

if (nargin < 2)
        error('Not enough input arguments.')
        return
end

for itype = 1:ntypes
  fprintf('Working on conversion of %s image...\n', img_types{itype});
  
  % Create the name of the first file in inDir
  firstfile = fullfile(inDir,sprintf('I.%03d', itype));

  % Generate the Analyze header and other info from the first file
  fprintf('Getting header info from: %s\n', firstfile);
  [header, orient, im_offset, adwcount] = GE_createSPMHeader(firstfile);
  volSize =  [header.dim(2) header.dim(3) header.dim(4)];
  
  % Reorient the header to SPM orientation
  header = GE_reorientHeader(header, orient); 

  imageVol = zeros(volSize);
  for islice = 1:volSize(3)
    % Determine name of current file
    curr_fname = fullfile(inDir,sprintf('I.%03d',(islice-1)*2+itype));
    
    % Open the file
    fid = fopen(curr_fname,'rb');
    if fid == -1
      fprintf('Failed to open: %s\n', curr_fname);
      imageVol = [];
      status = 0;
      return
    end
    
    % Read slice from file
    fseek(fid,im_offset+1,-1);
    buffer = fread(fid,volSize(1:2),sprintf('ubit%d',header.bitpix));
    
    % Append slice to image volume
    imageVol(:,:,islice) = reshape(buffer,volSize(1),volSize(2));
    
    % Close the file
    fclose(fid);
  end % for islice = 

  % Reorient the volume to SPM orientation
  imageVol = GE_reorientImage(imageVol, orient);

  outName = sprintf('%s_%s',outstub,img_types{itype});
  
  % Write the SPM header file
  outFile = strcat(outName,'.hdr');
  status = GE_writeSPMHeader(outFile,header);

  % Write out the Analyze File
  outFile = [outName sprintf('.img')];
  [fid,message] = fopen(outFile,'w');
  if (fid == -1)
    fprintf('Cannot Open %s for writing.\n',outFile);
    error(message);
  end

  fwrite(fid,reshape(imageVol,1,prod(size(imageVol))),'int16');
  status = fclose(fid);

  % Done with vol
  fprintf('Wrote %s.img.\n',outName);
end %for itype

% Done
fprintf('\nConversion Finished. \n');

return
