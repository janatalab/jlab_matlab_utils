function midi = read_mid(fname)
% midinfo = read_mid(fname);
%
% Read information in a text dump of a multitrack MIDI file
% The text dump is formatted by MIDI to Csound on a Mac.
%

% 3/18/01 PJ
%
% 8/3/01 PJ - changed from a structure with field arrays to a structure array
%             with scalar fields.  This will hopefully increase the speed of
%             reading .mid files.
%
% 08/09/01 PJ - Added fname to midi structure

if nargin < 1
  error('Please provide a filename')
end

fid = fopen(fname,'rt');
if fid < 0
  error(lasterr)
end

% Count the number of lines in the file and size the midi structure to
% accommodate all lines.  Delete unfilled lines later

disp('Counting the number of lines in the file');
num_lines = 0;
while ~feof(fid)
  fgetl(fid);
  num_lines = num_lines+1;
end
disp(sprintf('\nFound %d lines', num_lines))

% Initialize the midi structure and various arrays
disp('Creating structure to hold midi information ...')
%midi = struct('fname', '', 'track_names',{}, 'takt', [], 'onoff', [], ...
%    'pitch', [], 'vel', [], 'ctl_id', [], 'ctl_val', []);

track = zeros(num_lines,1);
takt = zeros(num_lines,1);
onoff = zeros(num_lines,1);
pitch = zeros(num_lines,1);
vel = zeros(num_lines,1);
ctl_id = zeros(num_lines,1);
ctl_val = zeros(num_lines,1);

disp(sprintf('Starting to read events'))
frewind(fid);
curr_track = 0;
ce = 0;
while ~feof(fid)
  lstr = fgetl(fid);
  [h count, errmsg, nextidx] = sscanf(lstr,'%s',1);
  switch(h)
    case 'Header:'
      targ_str = '#tracks';
      num_tracks = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
      disp(sprintf('Detected %d tracks', num_tracks))
      
      targ_str = 'division';
      qrt_gran = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
      disp(sprintf('Quarter note granularity: %d', qrt_gran))
    case '***'				% track start/stop marker
      is_start = findstr(lstr, 'start');
      if is_start
	curr_track = curr_track+1;
	disp(sprintf('Detected start of track %d', curr_track))
      end
    case 'Text'			% try to read track name
      lidx = findstr(lstr,'<') + 1;
      ridx = findstr(lstr,'>') - 1;
      track_names{curr_track} = lstr(lidx:ridx);
      disp(sprintf('Track name: %s', lstr(lidx:ridx)))
    otherwise
      if findstr(lstr, 'Tempo')
	tempo = sscanf(lstr(max(find(isspace(lstr))):end),'%d',1)/1000;
	disp(sprintf('Tempo (ms/qrt_note): %d', tempo))
      elseif findstr(lstr,'Note')
	ce = ce+1;

	track(ce) = curr_track;
	takt(ce) =  sscanf(lstr,'%d',1); % get timing
	
	if findstr(lstr,'Note on')
	  onoff(ce) = 1;
	else
	  onoff(ce) = 0;
	end
	
	targ_str = 'pitch';
	pitch(ce) = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
	targ_str = 'vel';
	vel(ce) = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
      elseif findstr(lstr, 'Ctl chg')
	ce = ce+1;
	track(ce) = curr_track;
	takt(ce) =  sscanf(lstr,'%d',1); % get timing
	
	targ_str = 'ctl';
	ctl_id(ce) = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
	targ_str = 'val';
	ctl_val(ce) = sscanf(lstr(findstr(lstr,targ_str)+length(targ_str):end),'%d',1);
	if ctl_val(ce)
	  onoff(ce) = 1;
	else
	  onoff(ce) = 0;
	end
      end
  end
end

disp(sprintf('Found a total of %d events', ce));

if ce < num_lines
  track(ce+1:end) = [];
  takt(ce+1:end) = [];
  onoff(ce+1:end) = [];
  pitch(ce+1:end) = [];
  vel(ce+1:end) = [];
  ctl_id(ce+1:end) = [];
  ctl_val(ce+1:end) = [];
end

disp('Creating final midi structure')
midi.fname = fname;
midi.tempo = tempo;
midi.qrt_gran = qrt_gran;
midi.track_names = track_names;
midi.track = track;
midi.takt = takt;
midi.onoff = onoff;
midi.pitch = pitch;
midi.vel = vel;
midi.ctl_id = ctl_id;
midi.ctl_val = ctl_val;

%midi.event_info = struct('track', track, 'takt', takt, 'onoff', onoff, ...
%    'pitch', pitch, 'vel', vel, 'ctl_id', ctl_id, 'ctl_val', ctl_val);

fclose(fid);
