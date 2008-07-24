function status = GE_convertADW(inDir,outStem,starttime,nvols)
%
% status = GE_convertADW(inDir,outStem,starttime,nvols)
%
%GE_convertADW
%
% Converts a series of GE slices acquired on the Advanced
% Development Workstation into Analyze format for SPM.
% inDir is the name of the first directory containing the
%      series, e.g. 003
% outStem is the stem used for naming the Analyze files
% starttime is the time point number to start
% nvols is the number of volumes (time points) to convert
%
% status = 1, error
% status = 0, all is well
%
% eg GE_convert('DATA/00023/003','converted',5,27)
% Will convert 27 time points starting at the 5th time point
% named converted_i0001.img, converted_i0002.img, ...,
% converted_i0027.img and their associated header files.
%
% Assumes the data is stored as 003/I.001 through 003/I.999 and
% then 023/I.001 through 023/I.999, etc.
% Remember to skip the template images ie 1 time point for the
% first run.  You don't need to do this for subsequent runs.
%
% Souheil J. Inati  
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%
% Modification history:
%
% 02/19/01 PJ Checks to see that volumes contain data, and quits the conversion
% if they don't

status = 0;

if (nargin < 4)
        error('Not enough input arguments.')
	status = 1;
        return
end

% Create the name of the first file in inDir
firstfile = fullfile(inDir,'I.001');

% Generate the Analyze header and other info from the first file
[header, orient, im_offset, adwcount] = GE_createSPMHeader(firstfile);
volSize =  [header.dim(2) header.dim(3) header.dim(4)];
% Reorient the header to SPM orientation
header = GE_reorientHeader(header, orient); 

% Loop over the volumes until there is no more looping to be done.
for i = 1:nvols
  passnum = starttime + i - 1;

  % Create analyze file name
%  vol_str = sprintf('000%d',i);
%  vol_str = vol_str(length(vol_str)-3:length(vol_str));
  outName = [outStem sprintf('_i%04d',i)];

  % Read in the volume
  [imageVol, lastfile] = GE_readVolume(inDir, passnum, volSize, ...
				       header.bitpix, im_offset);

  % Check to see that there is data in imageVol
  if isempty(imageVol)
    status = 1;
    break
  end
				   
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
    status = 1;
    error(message);
  end

  fwrite(fid,reshape(imageVol,1,prod(size(imageVol))),'int16');
  status = fclose(fid);

  % Done with vol
  if (status == 0)
    fprintf('Wrote %s\n',outName);
  else
    error('Error converting %s\n', outName);
  end

end

% Done
fprintf('\nConversion Finished. \n');

return
