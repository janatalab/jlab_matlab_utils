% extract_controller.m
%
% Extracts controller data

LOAD_DATA_FROM_MID = 0;
LOAD_DATA_FROM_TXT = 0;
LOAD_DATA_FROM_TR =0;

DO_PLOTS = 1;
PLOT_INDIV = 0;
PLOT_AVG = 1;

INITIAL_AMP_THRESHOLD = 5;
TIME_THRESHOLD = 2;

MIN_SEG_DUR = 0.1;  % minimum duration of 100 ms
MAX_SEG_DUR = 2;

MAX_GLUE_SEP = 0.5;  % max separation in seconds between segments to be glued

spline_tscale = 0:0.01:MAX_SEG_DUR;

% Columns in the pmat matrix
miditoolbox_params;

datapath = '/afs/cmb/usr/petr/projects/contpopout/data';
fstub = 'timbre_only_CC';
use_track = 2;

fname = fullfile(datapath,sprintf('%s.mid', fstub));

if LOAD_DATA_FROM_MID
  tr = javareadmidi(fname);
end

if LOAD_DATA_FROM_TXT
	tr = readmidi2tr(fname);
	outfname = fullfile(datapath, sprintf('%s_tr.mat',fstub));
	fprintf('Saving converted track info to %s\n', outfname);
	save(outfname,'tr')
end

if LOAD_DATA_FROM_TR
  infname = fullfile(datapath, sprintf('%s_tr.mat', fstub));
  fprintf('Loading data from %s\n', infname)
  load(infname)
end

% Parse individual events based on threshold

% Find all values above some threshold and find start and stop boundaries of
% each chunk 

time_vals = tr(use_track).pmat(:,CONTROLLER_SEC_COL);
slider_vals = tr(use_track).pmat(:,CONTROLLER_VAL_COL);

% Get potential segmentation markers based on time values
time_diff_vect = diff([0; time_vals]);
time_diff_mask = time_diff_vect > TIME_THRESHOLD;
time_seg_onset_idxs = find(time_diff_mask);
num_time_segs = length(time_seg_onset_idxs);

% Get potential segmentation markers based on controller values
diff_vect = diff([0; slider_vals>INITIAL_AMP_THRESHOLD; 0]);

% Check to see if the number of positive transitions matches the number of
% negative transitions
if sum(diff_vect)
  fprintf('Mismatched numbers of positive and negative transitions\n')
  seg_start_stop_idx = [];
else
  seg_start_stop_idx = [find(diff_vect>0) find(diff_vect<0)-1];
end

num_val_segs = size(seg_start_stop_idx,1);

% Display some information about segmentation markers that were identified
fprintf('Number of segments based on temporal threshold of %1.1f s: %d\n', TIME_THRESHOLD, num_time_segs);
fprintf('Number of segments based on value threshold of %d: %d\n', INITIAL_AMP_THRESHOLD, num_val_segs);

% Handle each individual segment that is above threshold
segment = struct('idxs',[],'tscale',[],'vals',[]);
for iseg = 1:num_val_segs
  start_idx = seg_start_stop_idx(iseg,1);
  stop_idx = seg_start_stop_idx(iseg,2);
  segments(iseg).startstop = [start_idx stop_idx];
  segments(iseg).tscale = time_vals(start_idx:stop_idx);
  segments(iseg).vals = slider_vals(start_idx:stop_idx);

  % make sure tscale values are unique
  bad_idxs = find(~diff(time_vals(start_idx:stop_idx)))+1;
  segments(iseg).tscale(bad_idxs) = [];
  segments(iseg).vals(bad_idxs) = [];
  
  segments(iseg).dur = max(segments(iseg).tscale)-min(segments(iseg).tscale);

  % some metrics for error-checking purposes
  segments(iseg).initial_is_max = segments(iseg).vals(1) == max(segments(iseg).vals);
  segments(iseg).final_is_max = segments(iseg).vals(end) == ...
      max(segments(iseg).vals);
  segments(iseg).nsamps = length(segments(iseg).vals);
  
end % for iseg=

% Flag segments as potentially bad based on various criteria
bad_dur_mask = [segments.dur] < MIN_SEG_DUR;
fprintf('\t%d segments shorter than minimum criterion of %1.2f s\n', sum(bad_dur_mask), MIN_SEG_DUR);

% Delete segments that are flagged as bad
segments(bad_dur_mask) = [];

% Identify that pairs of segments that should maybe be part of the same
% segment, i.e. they were severed based on sudden erroneous drop out that
% caused a boundary.  The following conditions are symptomatic of such cases:
%
% Bounds of segment are above some criterion threshold
%
% Segment severed during rise:
%   max(seg1) is final sample
%
% Segment severed during fall
%   max(seg2) is initial sample

max_final_mask = [segments.final_is_max];
max_initial_mask = [segments.initial_is_max];

glue_segs = find(all([max_final_mask(1:end-1); max_initial_mask(2:end)]));

fprintf('Found %d segments to glue together\n', length(glue_segs));
for iseg = 1:length(glue_segs)
  gs = glue_segs(iseg);
  
  % Make sure that segments to be glued are within a criterion time window of
  % each other
  seg_sep = segments(gs+1).tscale(1)-segments(gs).tscale(end);
  if seg_sep <= MAX_GLUE_SEP;
    fprintf('\tGluing segment %d and segment %d\n', [0 1]+gs);
    segments(gs).startstop(2) = ...
	segments(gs+1).startstop(2);
    segments(gs).tscale = [segments(gs).tscale; ...
	  segments(gs+1).tscale];
    segments(gs).vals = [segments(gs).vals; ...
	  segments(gs+1).vals];
    segments(gs).dur = diff(segments(gs).tscale([1 end]));
    segments(gs).nsamps = length(segments(gs).vals);
    segments(gs).final_is_max = 0;
  end
end % for iseg -- gluing

segments(glue_segs+1) = [];

% Toss segments that are too long
bad_dur_mask = [segments.dur] > MAX_SEG_DUR;
fprintf('\t%d segments longer than maximum criterion of %1.2f s\n', sum(bad_dur_mask), MAX_SEG_DUR);

% Delete segments that are flagged as bad
segments(bad_dur_mask) = [];

%
%  Interpolate the segments onto a common mesh so that they can be averaged
%

nseg = length(segments);
for iseg = 1:nseg
  % Set final value in segment to zero so that interpolation works
  segments(iseg).vals(end+1) = segments(iseg).vals(end);
  segments(iseg).tscale(end+1) = segments(iseg).tscale(end)+0.01;

  curr_tscale = segments(iseg).tscale-segments(iseg).tscale(1);
  
  % Do some spline interpolation and plot the waveform
  %segments(iseg).cs = spline(curr_tscale,[0; segments(iseg).vals; 0]);
  segments(iseg).cs = pchip(curr_tscale,[segments(iseg).vals]);
  
  segments(iseg).resamp = ppval(segments(iseg).cs,spline_tscale);
end


% Check for multiple peaks in a chunk. Multiple peaks might be an indication of
% multiple responses within a chunk, and/or the slider not being returned
% to zero.

if DO_PLOTS
  if PLOT_INDIV
    plot(curr_tscale,segments(iseg).vals,spline_tscale,segments(iseg).resamp,'r')
    title(sprintf('Segment %d/%d', iseg, nseg))
    pause
  end
  if PLOT_AVG
    figure(1),clf
    subplot(2,1,1)
    resp_mat = cat(1,segments.resamp);
    nresp = size(resp_mat,1);
    imagesc(spline_tscale,1:nresp,resp_mat), axis xy
    set(gca,'xlim',[0 MAX_SEG_DUR])
    
    subplot(2,1,2)
    errorbar(spline_tscale,mean(resp_mat),std(resp_mat)/sqrt(nresp))
    set(gca,'xlim',[0 MAX_SEG_DUR])
    title(sprintf('Mean response profile: %d responses',nresp))
  end
end