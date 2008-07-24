function GE_convert(inDir,outStem)
%
% GE_convert(inDir,outStem)
%
%GE_convert 
%
% Converts a series of GE slices into Analyze format for SPM.
% inDir is the name of the directory containing the series, e.g. 003
% outStem is the stem used for naming the Analyze files
%
% Souheil J. Inati  
% Dartmouth College 2000

%
% 04/18/00 PJ  Added check to make sure that imageVol is not empty when
% returning from GE_readVolume
%

if (nargin < 2)
        error('Not enough input arguments.')
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

% Loop over the runs
% Loop over the passes (volumes)
% Initialize the counters
done = 0;
runnum = 1;   % run number
passnum = 0;  % pas number

while ~done

  for i = 1:adwcount
    passnum = passnum + 1;

    % Create analyze file name
    series_str = sprintf('0%d',runnum);
    series_str = series_str(length(series_str)-1:length(series_str));
    vol_str = sprintf('000%d',i);
    vol_str = vol_str(length(vol_str)-3:length(vol_str));
    outName = [outStem sprintf('_s%s_i%s',series_str,vol_str)];

    % Read in the volume
    [imageVol, lastfile] = GE_readVolume(inDir, passnum, volSize, header.bitpix, im_offset);

    if isempty(imageVol)
      error('Failed to read in volume ...')
    end
    
    % Reorient the volume to SPM orientation
    imageVol = GE_reorientImage(imageVol, orient);

    % Write the SPM header file
    outFile = strcat(outName,'.hdr');
    status = GE_writeSPMHeader(outFile,header);

    % Write out the Analyze File
    outFile = [outName sprintf('.img')];
    [fid,message] = fopen(outFile,'w','l');
    if (fid == -1)
      fprintf('Cannot Open %s for writing.\n',outFile);
      error(message);
    end

    fwrite(fid,reshape(imageVol,1,prod(size(imageVol))),'int16');
    status = fclose(fid);

    % Done with vol
    fprintf('Wrote %s\n',outName);

    if runnum == 1 & adwcount > 1 & passnum == 1
       break
    end

  end % Run loop

  % Check for another run
  [path_stem, path_start] = fileparts(inDir);
  % Make the filename for the first file of the next run
  filep = ceil((passnum+1)*volSize(3)/999) - 1; 
  filen = (passnum+1)*volSize(3) - 999*filep;
  path_num = str2num(path_start) + 20*filep;
  path_now = sprintf('00%d',path_num);
  path_now = path_now(length(path_now)-2:length(path_now));
  path = fullfile(path_stem, path_now);
  stub = sprintf('00%d',filen);
  stub = stub(length(stub)-2:length(stub));
  nextFile = fullfile(path,['I.' stub]); % last file of next volume
  [fid,message] = fopen(nextFile,'r');

  if (fid == -1)
    done = 1;
  else
    fclose(fid);
    runnum = runnum + 1;
  end

end

% Done
fprintf('\nConversion Finished. \n');

return
