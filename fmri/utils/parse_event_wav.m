function [resp, events, slices] = parse_event_wav(fname,p)
%
% [resp, events, slices] = parse_event_wav(fname);
%
% Reads in a .wav file that contains stimulus and magnet event information.
% The default arrangement assumes that left and right button presses are on the
% left channel, and left button signals are smaller in magnitude than right
% button signals.  Magnet trigger and event onset information, together with
% the receiver unblank (slice timing) signals are mixed on the right.  Trigger
% and onset information are larger in magnitude than the receiver unblank.
%
% Button press onsets are negative going at start
%

% Originally written for prime_long

% 02/28/01 PJ Started adapting to be useful for other experiments, e.g. schub3
%             The function now only returns detected events on the two
%             channels, without any further elimination of events that might
%             occur on block boundaries.  That process is left to individual
%             scripts.
%
%             Implemented correct handling of button and event polarity

error(nargchk(1,2,nargin))
  
LOAD_DATA = 1;

if nargin < 2
  p = check_pinfo;				% check and fill pinfo
                                                % structure
else
  p = check_pinfo(p);
end

predicted_num_slices = p.slice_per_vol*(p.nvol+p.ndummy);

% Types of information on the different channels
info_type = {'button_press','stim_and_magnet'};

violation = [];
RESP = 1;
EVENT = 2;

LEFT = 1;
RIGHT = 2;

if LOAD_DATA
  disp(sprintf('\nLoading data from file: %s', fname))
  [data, Fs, nbits, opts] = wavread(fname);
end

nchan = size(data,2);

min_bp_samps = p.min_bp_isi/1000*Fs;
min_event_trig_samps = p.min_event_trig_isi/1000*Fs;

interslice_samps = p.interslice_intrvl/1000*Fs;

% Assume that the initial 1 sec contains no events, and use this to get
% estimate of mean and std of signal

basemean = mean(data(1:Fs,:));
basestd = std(data(1:Fs,:));

% Determine the maximum signal amplitude on each channel

if p.verbosity > 1
  disp(sprintf('Getting data min and max'))
end
[chanmax, chanmax_samps] = max(data);
[chanmin, chanmin_samps] = min(data);

if p.verbosity > 1
  for ichan = 1:size(data,2)
    disp(sprintf('Chan(%d): Minimum: %2.2f; Maximum: %2.2f', ichan,chanmin(ichan), chanmax(ichan)))
  end
end

%
%  HANDLE THE BUTTON PRESS CHANNEL
%

% Set thresholds for button presses
if p.button_polarity == -1
  thresh(1).high = 0.66*chanmin(1);	% right
  thresh(1).low = 0.25*chanmin(1);	% left
else
  if isfield(p,'button_threshhigh')
    thresh(1).high = p.button_threshhigh*chanmax(1);
  else
    thresh(1).high = 0.66*chanmax(1);	% right
  end
  if isfield(p,'button_threshlow')
    thresh(1).low = p.button_threshlow*chanmax(1);	% left
  else
    thresh(1).low = 0.25*chanmax(1);	% left
  end    
end

ichan = 1;

if p.verbosity > 1
  disp(sprintf('Processing response data on channel %d', ichan))
  disp(sprintf('Finding samples that exceed high (%1.2f) and low (%1.2f) thresholds', thresh(ichan).high, thresh(ichan).low))
end

switch info_type{ichan}
  case 'button_press'
    min_isi_samps = min_bp_samps;
  otherwise
    disp(sprintf('Do not know how to handle information type: %s', info_type{ichan}))
end

% Find all samples that exceed right and left button press threshold
if p.button_polarity == -1
  threshsamps_high = find(data(:,ichan) <= thresh(ichan).high);%right
  threshsamps_low = find(data(:,ichan) <= thresh(ichan).low); % left
else
  threshsamps_high = find(data(:,ichan) >= thresh(ichan).high);
  threshsamps_low = find(data(:,ichan) >= thresh(ichan).low); % left
end

% Get the first sample of each threshold crossing
cross_samps_high = threshsamps_high([1 (find(diff(threshsamps_high)>1)+1)']);

% Find all samps that exceed left button press threshold
cross_samps_low = threshsamps_low([1 (find(diff(threshsamps_low)>1)+1)']);

% Find all samps in low threshold vector that correspond to events in high
% threshold vector.  Allow several samples of slop in threshold determination

search_range = repmat(cross_samps_high',length(-p.perisamps:1:p.perisamps),1) + ...
    repmat((-p.perisamps:1:p.perisamps)',1,length(cross_samps_high));

search_range = search_range(:);
resp.samps = cross_samps_low;

resp.time = cross_samps_low/Fs;

rt_idx = find(ismember(cross_samps_low,search_range));
lft_idx = find(~ismember(cross_samps_low,search_range));
resp.id(rt_idx) = RIGHT;
resp.id(lft_idx) = LEFT;

disp(sprintf('Found a total of %d threshold crossings: RIGHT: %d; LEFT: %d', length(resp.samps), length(rt_idx), length(lft_idx)))

% Make sure that no putative events occur within the minimum button press ISI
tmpdiff = diff(resp.samps);
bad_idx = find(tmpdiff < min_bp_samps)+1;
num_bad = length(bad_idx);  

% Eliminate overlapping button press events
if num_bad
  if p.verbosity
    disp(sprintf('Found %d overlapping button press events', num_bad))
    disp('Eliminating ...')
  end
  resp.samps(bad_idx) = [];
  resp.time(bad_idx) = [];
  resp.id(bad_idx) = [];
end

% Do some sanity checks
if length(resp.samps) ~= p.expect_num_resp
  if p.verbosity
    disp(sprintf('     WARNING: Number of detected RESPONSES (%d) did not match expected number (%d)', length(resp.samps), p.expect_num_resp))  
  end
  violation(RESP) = 1;
end

if nargout == 1
  return
end

%
%  HANDLE THE EVENT ONSET AND MAGNET TRIGGER CHANNEL
%

ichan = 2;

if p.verbosity > 1
  disp(sprintf('\nProcessing stimulus event data on channel %d', ichan))
end

%
%  Deal with the event signals first
%

min_isi_samps = min_event_trig_samps;

% Set threshold for triggers

if p.event_thresh & p.event_polarity
  thresh(ichan).high = p.event_thresh * chanmax(ichan) * p.event_polarity;
else
  thresh(ichan).high = 0.5;
end

if p.verbosity > 1
  disp(sprintf('Finding samples that exceed high (%1.2f) threshold', thresh(ichan).high))
end

if p.event_polarity == -1
  threshsamps_high = find(data(:,ichan) <= thresh(ichan).high);
else
  threshsamps_high = find(data(:,ichan) >= thresh(ichan).high);
end

cross_samps_high = threshsamps_high([1 (find(diff(threshsamps_high)>1)+1)']);
bad_idx = find(diff(cross_samps_high(:)) < min_isi_samps);

disp(sprintf('Found a total of %d threshold crossings', length(cross_samps_high)));

num_bad = length(bad_idx);

if num_bad
  disp(sprintf('Found %d bad event triggers.  Eliminating ...', num_bad))
  cross_samps_high(bad_idx) = [];
end

%
% Realign everything to first event (assuming this is the magnet trigger
%

if ~isfield(p,'assume_first_is_magtrig')
  p.assume_first_is_magtrig = 1;
end

if p.assume_first_is_magtrig
  if p.verbosity > 1
    disp(['   Stripping first event (run trigger) and calculating event times' ...
	  ' relative to it ...'])
  end

  run_startsamp = cross_samps_high(1);
  event_samps = cross_samps_high(2:end)-run_startsamp;

  % Remove sample offset due to dummy shots
  ndummy_samps = p.ndummy*2*Fs;
  event_samps = event_samps - ndummy_samps;
  events.samps = event_samps;
  events.time = events.samps/Fs;

  % Subtract start offset from reaction time data
  resp.samps = resp.samps-run_startsamp - ndummy_samps;
  resp.time = resp.samps/Fs;

else
  events.samps = cross_samps_high;
  events.time = events.samps/Fs;
end

return

%
% DO NOT DO ANY FURTHER PROCESSING OF EVENT DATA.  NOW THAT EVENTS HAVE BEEN
% IDENTIFIED, LET CALLING SCRIPT TO FURTHER PROCESSING.  
%


% Eliminate the events corresponding to noise bursts and beginning and end of
% block boundaries

if p.verbosity > 1
  disp('   Eliminating noise burst markers ...')
end
tmp = find(diff(event_samps) > 20*Fs)+1;
tmp = [tmp tmp-1];
tmp = tmp(:);

if max(tmp) ~= length(event_samps)
  event_samps([1 tmp' length(event_samps)]) = [];
else
  event_samps([1 tmp']) = [];
end

events.samps = event_samps;
events.time = events.samps/Fs;

% Check to see if number of event samps matches the expected number
if length(event_samps) ~= p.expect_num_events
  if p.verbosity
    disp(sprintf('     WARNING: Number of detected EVENTS (%d) did not match expected number (%d)', length(event_samps), p.expect_num_events))
  end
  violation(EVENT) = 1;
end

if nargout < 3
  return
end

%
%  Deal with the slice timing info
% 

thresh(ichan).low = 0.2 * chanmin(ichan); % threshold for slices

threshsamps_low = find(data(:,ichan) <= thresh(ichan).low);
cross_samps_low = threshsamps_low([1 (find(diff(threshsamps_low)>1)+1)']);

% Get rid of first threshold crossing as this is the off edge of the start
% trigger

cross_samps_low(1) = [];

% Toss out any negative crossings that occur earlier than known slice period

tmp = [cross_samps_low cross_samps_low+interslice_samps-p.interslice_slop_samps]';
tmpdiff = diff(tmp(:));

badidx = find(tmpdiff(2:2:end) < 0)+1;
longidx = find(tmpdiff(2:2:end) > (2*p.interslice_slop_samps)) ...
    + 1;

% Some of the samples in badidx are actually good.  This can be determined by
% the proximity of indices given in badidx.  If they are adjacent samples,
% chances are very good that the first one represents the off-edge of an event
% trigger and the following one represents the onset of the next slice.
% Therefore, if two indices differ by only 1, then toss the first.
% Sometimes the off-edge of an event trigger masks the onset of the next slice,
% and the delay until the next slice onset is longer than the period between
% slice onsets. In these cases, slice onsets have to be inferred and inserted
% into the sequence.

% Handle the case of missed slice onsets first.  The simplest fix is to replace
% the sample numbers of the extra short intervals with corrected indices. 

if ~all(ismember(longidx-1,badidx))
  if p.verbosity
    disp('WARNING: Extra long duration not preceded by extra short duration')
  end
end

cross_samps_low(longidx-1) = mean(cross_samps_low([longidx-2 longidx]),2);
badidx(find(ismember(badidx,longidx-1))) = [];

badidx(find(diff(badidx)==1)+1) = [];

cross_samps_low(badidx) = [];		% remove bad samples

% Reiterate process to toss out any negative crossings that occur earlier than known slice period
tmp = [cross_samps_low cross_samps_low+interslice_samps-p.interslice_slop_samps]';
tmpdiff = diff(tmp(:));

badidx = find(tmpdiff(2:2:end) < 0)+1;
badidx(find(diff(badidx)==1)+1) = [];
cross_samps_low(badidx) = [];		% remove bad samples

if predicted_num_slices ~= length(cross_samps_low)
  error(sprintf('Predicted number of slices (%d) does not match number found (%d)', predicted_num_slices, length(cross_samps_low)))
end

function p = check_pinfo(p)
  if nargin < 1
    p = struct(...
	'verbosity', [], ...
	'min_bp_isi',[], ...
	'min_event_trig_isi',[], ...
	'perisamps', [], ...
	'interslice_slop_samps', [], ...
	'interslice_intrvl', [], ...
	'slice_per_vol', [], ...
	'nvol', [], ...
	'ndummy', [], ...
	'event_polarity', [], ...
	'event_thresh', [], ...
	'expect_num_events', [], ...
	'slice_polarity', [], ...
	'button_polarity' ...
	);
  end
  
  if isempty(p.verbosity)
    p.verbosity = 1;
  end
  
  if isempty(p.min_bp_isi)
    p.min_bp_isi = 2000;			% minimum duration between button presses (ms)
  end
  
  if isempty(p.min_event_trig_isi)
    p.min_event_trig_isi = 1000;
  end
  
  if isempty(p.perisamps)
    p.perisamps = 5;
  end
  
  if isempty(p.interslice_slop_samps)
    p.interslice_slop_samps = 10;
  end

  if isempty(p.interslice_intrvl)
    p.interslice_intrvl = 74;			% time between slice onsets (ms)
  end
					
  if isempty(p.slice_per_vol)
    p.slice_per_vol = 27;
  end
  
  if isempty(p.nvol)
    p.nvol = 184;
  end
  
  if isempty(p.ndummy)
    p.ndummy = 4;
  end
  
  if isempty(p.event_polarity)
    p.event_polarity = 1;
  end
  
  if isempty(p.event_thresh)
    p.event_thresh = 0.4;
  end
  
  if isempty(p.expect_num_events)
    p.event_thresh = 0;
  end
  
  if isempty(p.slice_polarity)
    p.slice_polarity = -1;
  end
  
  if isempty(p.button_polarity)
    p.button_polarity = -1;
  end
  
 % end check_pinfo
