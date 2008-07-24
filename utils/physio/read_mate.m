function [ri, p] = read_mate(flist,p)
% [runinfo, params] = read_mate(flist,params)
%
% Reads physiological data monitoring files acquired using the Siemens 3T MATE
% system.
%
% Each variable is recorded to a different file, and these files are specified
% in flist.  There is a field in the params structure which identifies the filetypes

ri = struct('vol_onsets',[],'cardiac',[]);

if nargin < 2
  error('read_mate: Too few input arguments')
end

% Determine how many files we are dealing with
nfiles = length(flist);

% Make sure we have a field specifying their types
if ~isfield(p,'ftypes') || length(p.ftypes) ~= nfiles
  fprintf('read_mate: File type information incorrectly specified');
end

% One of the files must be the external trigger signal from the scanner, as
% this will be used to determine start and stop times of a run, and for parsing
% of the other files

% Load in all of the files
clear trig pulse resp
for ifile = 1:nfiles
  fname = flist{ifile};
  
  if ~exist(fname,'file')
    fprintf('read_mate: Could not find file <%s> ... skipping ...\n', fname);
    continue
  end
  
  fprintf('Loading data from file: %s\n', fname);
  data = load(fname, '-ascii');
  
  switch p.ftypes{ifile}
    case {'scantrig', 'exttrig'}
      trig = data;
    case {'cardiac', 'pulse'}
      pulse = data;
    case {'respir'}
      resp = data;
    otherwise
      fprintf('read_mate: Unknown filetype: %s\n', p.ftypes{ifile});
  end
end % for ifile=
clear data

if exist('trig')
  % Process the trigger data
  tp.thresh = p.trig_thresh;
  tp.dir = p.trig_dir;

  % Find the trigger pulse onsets
  trig_onset_samps = find_thresh_cross(trig,tp);

  % See how many trigger onsets we found, and make sure this number is reasonable
  ntrig = length(trig_onset_samps);
  if ~ntrig
    fprintf('read_mate: Found no scanner triggers ...\n')
    return
  end

  % Convert the onsets to time values
  trig_onset_sec = (trig_onset_samps-1)/p.Fs_trig;

  % Calculate some statistics about this pulse train that will be used to
  % determine run onsets and offsets.

  % Get the median inter-trigger interval and make sure it is within bounds
  itis = diff(trig_onset_sec);
  median_iti = median(itis);

  if p.trig_is_vol
    expect_median_iti = p.TR;
  else
    expect_median_iti = p.TR/p.nslice_per_vol;
  end
  min_iti = expect_median_iti-p.trig_slop_sec;
  max_iti = expect_median_iti+p.trig_slop_sec;
  if ~((median_iti > min_iti) & (median_iti < max_iti))
    warning(sprintf('read_mate: median iti (%1.4f s) is out of range (%1.4f, %1.4f)\n', median_iti, min_iti, max_iti))
    return
  end

  % Find the ITI that exceeds some integer multiple of the median ITI
  border_idxs = find(itis > 3*median_iti)+1;
  num_borders = length(border_idxs);

  % If we find none, then assume that all of the trigger events are part of the
  % same run
  if num_borders == 0
    trig_onset_idxs = 1:ntrig;
  elseif num_borders == 1;
    trig_onset_idxs = border_idxs:ntrig;
  else
    fprintf('Located %d run borders.  Too many to handle\n', num_borders);
    return
  end

  run_start_time = trig_onset_sec(trig_onset_idxs(1));
  run_stop_time = trig_onset_sec(trig_onset_idxs(end))+p.TR;
  ri.vol_onsets = trig_onset_sec(trig_onset_idxs)-run_start_time;
else
  fprintf('Have no trigger information: Cannot process other information\n');
  return
end % if exist('trig')
  
for ifile = 1:nfiles
  switch p.ftypes{ifile}
    case {'cardiac', 'pulse'}
      tp.thresh = p.pulse_thresh;
      tp.dir = p.pulse_dir;
      
      % Get the pulse signal threshold crossings
      pulse_onset_samps = find_thresh_cross(pulse,tp);
      pulse_onset_sec = pulse_onset_samps/p.Fs_pulse;
      
      % Get the pulse events for the run
      pulse_idxs = find((pulse_onset_sec >= run_start_time) & (pulse_onset_sec ...
	  <= run_stop_time));
      
      ri.cardiac = pulse_onset_sec(pulse_idxs)-run_start_time;
      
  end % switch p.ftypes{ifile}
end % for ifile

return

