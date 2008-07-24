function status = GE_convertVolume(inDir,runnum,outName)
%
% status = GE_convertVolume(inDir,runnum,outName)
%
%GE_convertVolume
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
% 08/30/04  Petr Janata -- added some filename handling (checking for alternate firstfile)

status = 1;

if (nargin < 3)
        error('Not enough input arguments.')
        return
end

% Create the name of the first file in inDir
firstfile = fullfile(inDir,'I.001');

% On rare occassions, it might be the case that the I is lowercase.  Check to
% see if this is the case
if ~exist(firstfile)
  new_first = fullfile(inDir,'i.001');
  if exist(new_first)
    fprintf('Expected I.001, found i.001.  Using i.001\n');
    firstfile = new_first;
  else
    error(sprintf('Could not locate first file in image: %s\n', firstfile));
  end
end 

% Generate the Analyze header and other info from the first file
[header, orient, im_offset, adwcount] = GE_createSPMHeader(firstfile);
volSize =  [header.dim(2) header.dim(3) header.dim(4)];
% Reorient the header to SPM orientation
header = GE_reorientHeader(header, orient); 

% Read in the volume
[imageVol, lastfile] = GE_readVolume(inDir, runnum, volSize, header.bitpix, im_offset);

% Reorient the volume to SPM orientation
imageVol = GE_reorientImage(imageVol, orient);

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

% Done
fprintf('\nConversion Finished. \n');

return
