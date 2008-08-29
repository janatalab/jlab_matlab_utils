function [channels, bfactors, npts, dt] = parseheader(topdir, hdrfile)
% PARSEHEADER obtain study information from Bti header file
%
% Read header file (output of Bti 'convert') for information necessary to
% read Bti data files directly without requiring the database.

if ~exist('hdrfile') hdrfile = []; end
if isempty(hdrfile)			% determine header file name
  dd = dir(topdir);
  for i=1:length(dd)
    if findstr([dd(i).name 'zzz'],'hdr')
      hdrfile = dd(i).name
    end
  end
end
if isempty(hdrfile)
  error('Cant find header file');
end
fid = fopen([topdir hdrfile],'rt');	% Open header file

line = 0;
while line ~= -1			% parse ascii lines within file
  tmpline = [line 'zzz'];
  if findstr(tmpline, 'Total Channels')
    nchans = str2num(line((findstr(line,': ')+2):end));
  elseif findstr(tmpline, 'Sample Period')
    dt = str2num(line((findstr(line,': ')+2):end));
  elseif findstr(tmpline, 'Total Points')
    npts = str2num(line((findstr(line,': ')+2):end));
  elseif findstr(tmpline, 'MxA')
    channels = line;
  elseif findstr(tmpline, 'MSI.Conversion')
    bfactors = line;
  end
  line = fgets(fid);
end

fclose(fid);

rchans = channels; rfacts = bfactors;
channels = []; factnames = [];
[junk rfacts] = strtok(rfacts);		% readout line label first
for i=1:nchans
  [tchan rchans] = strtok(rchans);	% parse channel names
  channels = strvcat(channels, tchan);
  [tfact rfacts] = strtok(rfacts);	% parse converstion factor strings
  factnames = strvcat(factnames, tfact);
end

bfactors = str2num(factnames);
