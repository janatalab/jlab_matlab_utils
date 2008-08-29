function [B, channel, npts, dt] = getbtiraw(topdir,datfile,hdrfile)
% GETBTIRAW get MEG data from Bti raw data files
%
% Open file topdir/datfile, using hdrfile to translate the rawfile format.
% Return magnetic field data B (zero-meaned), channel names (both suitably
% reordered), number of points npts, and time step dt.
% At present doesn't return n_epochs or initial time.

% Example: datfile = '3047.4/2Mx3-6s254/GV8RH9BO/1/c,rfhp1.0Hz'

% Get format information from header file
if ~exist('hdrfile') hdrfile = []; end
[channel bfactors npts dt] = parseheader(topdir, hdrfile);

fid = fopen(fullfile(topdir,datfile), 'r', 'ieee-be');
% use integer*2 format if reading raw digitized data (short int)
if ~exist('rformat'), rformat = 'integer*2'; end

% else data may be in float format if filtered or DC offset

nchans = size(channel,1);
[B num] = fread(fid,[nchans npts],rformat);
if num ~= (npts*nchans)
  error('Number of data points read doesnt match number requested')
end

fclose(fid);

sortarray = [];				% sort channels into good order
for i=1:nchans
  if channel(i,1) == 'A'		% good MEG channels first
    sortarray(i) = str2num(channel(i,2:end));
  elseif channel(i,1:3) == 'TRI'	% Trigger next
    sortarray(i) = 500;
  elseif channel(i,1:3) == 'RES'	% then Response channel
    sortarray(i) = 501;
  elseif channel(i,1) == 'X'		% finally external channels
    sortarray(i) = 1000 + str2num(channel(i,2:end));
  else
    sortarray(i) = 2000 + i;
  end
end

[channum sortindex] = sort(sortarray);
sortindex = sortindex(channum<2000);	% discard extraneous channels
B = B(sortindex,:)';
channel = channel(sortindex,:);
bfactors = bfactors(sortindex);
bmeans = mean(B);
nbchans = find(strcmp(cellstr(channel),'TRIGGER')) - 1; % # of Bfield channels
bmeans(nbchans+1:end) = 0*bmeans(nbchans+1:end); % only zeromean B channels

if rformat == 'integer*2'		% apply channel factors to ints
  for i=1:size(B,2)
    B(:,i) = bfactors(i)*(B(:,i)-bmeans(i));
  end
end

