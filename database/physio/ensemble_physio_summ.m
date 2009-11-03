function outdata = ensemble_physio_summ(indata,defs)

% calculates summary statistics for physio data
% 
% This function calculates a variety of summary statistics for epochs of
% physiological data signals, matches those statistics with responses from
% presentation logfiles (optional), saves these data to disk (signals and
% peaks) in a form amenable to import into SAS (optional), and also saves
% waveform graphics that may be used in assessing the validity of the
% summary statistics that have been calculated, as well as the processing
% steps taken during those calculations.
% 
% This script has the option of either completing everything in one step,
% or using a two-step process whereby epochs and peaks are extracted in the
% first step, and then a second step takes the epoched and peaked data and
% calculates summary statistics from them. This allows manual intervention
% between steps to visually inspect and remove bad peaks or waveform
% artifacts from the data.
% 
% NOTE: this script assumes that pre-processed physio data have been saved
% into eeglab .set files, and that presentation log files are available to
% guide signal epoching
% 
% REQUIRES
% 
%   indata - cell array of ensemble data structs
%       'sinfo' - data struct, containing sinfo data for all (but only)
%           those subjects that you wish to analyze
%       'paths' - ensemble data struct, output of ensemble_fmri_init,
%           containing 'physio_outdir', which is used as the output
%           location of all figures that are constructed
%       'physio_paths' - output of ensemble_physio, ensemble data
%           struct containing one row per subject/session/run, with paths
%           to pre-processed physio data files for each run
%       'pres_paths' - output of ensemble_parse_present, an ensemble
%           data struct containing one row per subject/session/run, with
%           paths to pre-processed presentation response and event data
% 
% LOADING PHYSIO DATA FILES
%   
%   physio files are loaded on a per-subject, per-session, per-run basis,
%   from locations identified in 'physio_paths', for those subjects that
%   exist in 'sinfo'. By default, ensemble_physio (which generates
%   physio_paths) generates a non-filtered raw data file for each
%   subject/session/run, and saves these to physio_paths, and given
%   filtering criteria, also saves a filtered data file to disk in the same
%   location as the raw data file, and saves the path to the filtered data
%   file into physio_paths, with '_filtered' appended to the end of the
%   filename, before the extension. If the following parameter is set to 1,
%   ensemble_physio_summ will search for and attempt to load the
%   filtered data file as opposed to the raw data file. If set to 0, it
%   will attempt to load the raw, unfiltered data file.
% 
%   defs.physio_summ.use_filtered (default 1)
% 
% EPOCHING
% 
%   ensemble_physio_summ is intended to calculate summary statistics
%   on epochs of physiological data signals. These epochs may, or may not,
%   have responses or experimental classes associated with them. Epochs can
%   currently be defined in terms of events from a presentation log file,
%   or vectors of epoch onsets and offsets in seconds (calculated from the
%   beginning of the provided preprocessed physiological signal)
% 
%   defs.physio_summ.epoch_start (default 0) - if a struct, with sub-field
%       'filter', this is treated as a filter criteria for the included
%       presentation data, to filter for events to be treated as epoch
%       onsets. If a numeric vector, treated as a vector of epochs onsets
%       in seconds. If zero, no epoching will be conducted.
%   defs.physio_summ.epoch_end (default 0) - if single-length numeric, this
%       is treated as the epoch length in seconds. If a numeric vector of
%       the same length as epoch_start or the onset vector returned by
%       epoch_start, this is treated as a vector of epoch offsets in
%       seconds. If a struct, with sub-field 'filter', this is treated as
%       filter criteria for the included presentation data, to filter for
%       events to be treated as epoch offsets. By default, epoch_starts
%       will be matched with epoch_ends, and epoching will only continue if
%       an equal number of starts and ends is found, with all matching ends
%       occuring after their supposed matched start. By default, the max
%       epoch length will be used if a vector is found as the result of a
%       filter. If a filter criterion is used and the field 'minimum_end'
%       exists and is set to = 1, the minimum epoch time (offset-onset)
%       will be used for all epochs. If this variable is set to 0 and
%       defs.physio_summ.epoch_start is not set to 0, epoch lengths will be
%       set to the maximum interval between epoch starts. NOTE: this will
%       cause overlapping signal between some epochs where there is a
%       variance in the inter-start interval. If this is a problem, you
%       should set 'minimum_end' = 1, so that the epoch length gets set to
%       the minimum inter-start interval.
%   defs.physio_summ.baseline (default 2) - pre-epoch baseline to be added to
%       epoch when calculating and extracting epochs. this value is assumed
%       to be in 'seconds'.
% 
%       FIXME: In the case that there is not sufficient room
%       for the given baseline period between the first epoch onset and the
%       very beginning of the signal acquisition (for instance, first epoch
%       onset is at 1s, but the baseline period is 2s), the first epoch
%       will be not be returned by eeglab. This is because eeglab
%       requires an equal number of datapoints for all epochs extracted
%       from a given signal, and this case creates an out-of-bounds epoch
%       definition that is dropped by eeglab.
% 
%       FIXME: pad the signal with 0s to provide enough data points to
%       allow for a baseline for the first epoch, and then adjust timings
%       of all other events? this would necessarily require special
%       calculation of the baseline of the first epoch, limited to
%       real-data only, when subtracting the baseline from the signal of
%       interest, so as not to bias the baseline with artifically added 0s.
% 
% SUMMARY DATA CALCULATION
% 
%   The current analyses calculate summary statistics for both pulse and
%   skin conductance data. For pulse, mean heart rate within an epoch, as
%   well as heart rate variance, is calculated. For scr, a number of
%   statistics are calculated to investigate both individual transient skin
%   conductance responses (SCRs), and change in overall skin conductance
%   level (SCL).
% 
%   For SCRs, individual responses are defined as the highest peak within a
%   range of +/- 'proxlim' seconds (default 0.6) within an epoch, whose
%   height is greater than height_thresh (default 0.02). A standard minimum
%   SCR height, as defined by Dawson et al 2000, is 0.02 microsiemens,
%   therefore if the incoming signal is defined in microsiemens (as is
%   biopac data when both high pass filters on the GSR100 are set to DC),
%   the default height_thresh should be used. Total number of SCRs within
%   an epoch (npeaks), mean peak height, peak height variance, first peak
%   height, and first peak time are calculated. If preprocessed data are
%   provided in microsiemens, then heights and variance will be returned in
%   microsiemens.
% 
%   For SCLs, epoched skin conductance data are standardized
%   within-subject, across epochs, to control for individual differences in
%   SCL variability and range (Dawson et al 2000, p209; also Ben-Shakhar
%   1985). Specific integral and mean SCL over the expected period of a
%   specific SCR (spc_integral.being, spc_integral.end, default 2s to 6s)
%   are calculated, as well as integral and mean SCL over the entire epoch
%   (not including the baseline period).
% 
%   NOTE: if entering these data into mixed effects multi-level regression
%   models, certain effects such as different rates of habituation between
%   individuals will be somewhat controlled for by pooling variance within
%   subjects, but if not, additional measures should be considered if you
%   want to control for individual differences in (at least) rate of
%   habituation when using repeated stimuli or chronic stimuli. Different
%   rates of habituation may effect mean peak height and peak height
%   variance, and may also have some effect on integral.
% 
%   defs.physio_summ.channels - channel names, as they appear in your
%       eeglab chanlocs labels, that you wish to calculate summary
%       information for.
%   defs.physio_summ.scr - settings for SCR/GSR data summary
%       .spc_integral.begin (default 2) - time, in seconds, to start
%           calculation of the specific SCR integral
%       .spc_integral.end (default 6) - time, in seconds, to end
%           calculation of the specific SCR integral
%       .find_peaks - the following parameters are sent to find_peaks:
%   		.thresh (default 0) - proportion of total signal
%               range below which peaks will not be considered. for
%               example, for .thresh = 0.5, peaks will only be returned if
%               their height is more than half of the total signal range.
%               This is passed on to find_peaks, but is always set to 0,
%               since we want to be able to find peaks across the entire
%               range of signal.
%           .peak_heights.calcFrom (default 'left') - which side of a peak
%               to find the trough used to calculate peak height
%           .peak_heights.transformation (default 'none') - the
%               transformation, if any, performed on peak height values
%               before returning them.
%           .peak_height_window (default 0.6) - sliding window in
%               which only the highest peak is retained. for each peak
%               found, if there are other peaks within +/- seconds of that
%               peak, all but the highest peak will be discarded.
%   		.peak_height_thresh (default 0.02) - peak heights below
%               this value will be rejected
%       .keep_baseline (default: 1) - keep the baseline period in the
%           final signal that is returned (1), or remove it (0)? either
%           way, the baseline period (if greater than 0) mean is removed
%           from the signal of interest for each epoch, and peaks during
%           the baseline period are removed after peak finding
%   defs.physio_summ.pulse - settings for pulse ox. data summary
%   	.find_peaks.thresh (default 0.5) - proportion of total signal
%           range below which peaks will not be considered. for example,
%           for .thresh = 0.5, peaks will only be returned if their height 
%           is more than half of the total signal range
% 
% RESPONSE MATCHING FOR EPOCHS
% 
%   if there are responses expected and present within a presentation
%   logfile that can be matched in a 1 to 1 relationship with epochs that
%   are extracted or defined, these responses can be included in any data
%   files that are exported, matched to each epoch and the summary
%   statistics calculated for each epoch. Also, overlays of all epochs
%   corresponding to each unique response within a response set, and
%   average waveform for each level of each response, can be calculated and
%   output. presentation data must be included, in the form of a pres_paths
%   data struct (created by presentation_parse) containing paths to
%   presentation data that has been parsed and saved to disk for the given
%   participant/session.
% 
%   defs.physio_summ.responses (default 0) - struct array, each struct
%       identifying filtering criteria to be applied to the given
%       presentation data, to extract events containing response data that
%       should be mapped to each trial. This response data will be exported
%       along with summary information calculated for each trial.
%       .type - string, name of the given response. this string will be
%           used as the variable name for the given response data in the
%           resulting output data struct
%       .event - struct, assumed to include either an .include or .exclude
%           fieldname, and subsequent ensemble_filter criteria for the
%           given response class
% 
% DATA EXTRACTION VS SUMMARY STATISTIC CALCULATION
% 
%   the default is to extract epochs, extract peaks, and calculate summary
%   statistics in one step. You can control the execution of each step with
%   the following parameters:
% 
%   defs.physio_summ.EXTRACT_DATA (default: 1) - extract epochs and peaks
%       from provided waveforms. These data are returned in 'gsr_epochs'
%       and 'cardiac_epochs' datasets. REQUIRES: physio_paths and
%       pres_paths input datasets
%   defs.physio_summ.CALCULATE_STATS (default: 1) - calculate summary
%       statistics on extracted epochs and peaks. REQUIRES: either
%       defs.physio_summ.EXTRACT_DATA, physio_paths, and pres_paths, OR
%       'gsr_epochs' and 'cardiac_epochs' input data structures
% 
%   to use the two step process, you can first execute this script with
%   EXTRACT_DATA=1, CALCULATE_STATS=0, then import the results from that
%   step into another script (possibly manual_adjust_signal_data.m) to
%   manually edit the waveforms and peak information, and then return those
%   edited datasets to this script with EXTRACT_DATA=0, CALCULATE_STATS=1
%   to complete the summary statistic calculation.
% 
% DATA OUTPUT
% 
%   defs.physio_summ.GRAPH_EPOCHS (default 0)
%   defs.physio_summ.EPOCHS_BY_RESP (default 0)
%   defs.physio_summ.SAVE_EPOCHS (default 1)
%   defs.physio_summ.SAVE_DATA (default 0)
%   defs.physio_summ.export (default 0)
%   defs.physio_summ.sas (default 0)
%   defs.physio_summ.export_stub (default '')
% 
% FIXME - should be re-written, to include a mechanism by which channels to
% analyze are identified by EEG.chanlocs name as sub-fields of the params
% struct, and under that, sub-fields of each channel should be named
% exactly as the names of functions used to either filter or report the
% data. Then, sub-fields of those fields should be params that get sent to
% those functions as the second argument, the first argument being the
% data output of the previous function (in may or all cases, an eeglab
% struct? ... OR a cell array of structs containing a series of processing
% jobs, params for those jobs, function handles/references, etc ...?
% 
% CITATIONS
%   Dawson ME, Schell AM, & Filion DL (2000). "The Electrodermal System."
%       In JT Cacioppo, LG Tassinary & GG Berntson (Eds.) "Handbook of
%       Psychophysiology" (pp 201-223). Cambridge: Cambridge Univ. Press.
%   Ben-Shakhar  (1985). "Standardization within individuals: Simple method
%       to neutralize individual differences in skin conductance."
%       Psychophysiology, 22:292-9.
% 
% FB 2009.04.21
% FB 2009.05.04 - adding documentation, now zscores signal for SCL
% calculation but leaves signal for SCR calculation in raw form.
% FB 2009.06.18 - changed "params.lead" to "params.baseline", and checked
% the script for hard-coded reliance on anything fmri-related, with aims to
% generalize from ensemble_fmri_physio_summ to ensemble_physio_summ

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'physio_summ';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% % load colormap
% load new_seismic;
% colormap(new_seismic);

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
      case {'physio_paths'}
        physio_paths = indata{idata};
        physcol = set_var_col_const(physio_paths.vars);
      case {'pres_paths'}
        pres_paths = indata{idata};
        prescol = set_var_col_const(pres_paths.vars);
      case {'sinfo'}
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case {'gsr_epochs'}
        gsr_epochs = indata{idata};
        gecol = set_var_col_const(gsr_epochs.vars);
      case {'cardiac_epochs'}
        cardiac_epochs = indata{idata};
        cecol = set_var_col_const(cardiac_epochs.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','pathdata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && ~isempty(pathdata.data{1})
    if exist('nsub_proc','var') && length(nsub_proc) == 1 && ...
            exist('proc_subs','var')
      pfilt = struct();
      pfilt.include.all.subject_id = proc_subs;
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        sfilt = pfilt;
        sfilt.include.all.path_type = {'physio_outdir'};
        spathdata = ensemble_filter(lpathdata,sfilt);
        if length(spathdata.data{1}) == 1
          % one epi outdir, save outdata = epi_outdir
          outdata = spathdata.data{pcol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pcol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pcol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% initialize outdata struct for sinfo
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% check for physio_summ parameters within defs
if isfield(defs,'physio_summ')
  ps = defs.physio_summ;
else
  error('please provide params through defs.physio_summ\n');
end

% get parameters
try chans_req = ps.channels; catch error('must provide channel names\n'); end
try GRAPH_EPOCHS = ps.GRAPH_EPOCHS; catch GRAPH_EPOCHS = 0; end
try EPOCH_BY_RESP = ps.EPOCH_BY_RESP; catch EPOCH_BY_RESP = 0; end
try SAVE_EPOCHS = ps.SAVE_EPOCHS; catch SAVE_EPOCHS = 0; end
try SAVE_DATA = ps.SAVE_DATA; catch SAVE_DATA = 0; end
try EXTRACT_DATA = ps.EXTRACT_DATA; catch EXTRACT_DATA = 1; end
try CALCULATE_STATS = ps.CALCULATE_STATS; catch CALCULATE_STATS = 1; end
try export = ps.export; catch export.write2file = 0; end
try sas = ps.sas; catch sas.libsave = 0; end
try epoch_start = ps.epoch_start; catch epoch_start = ''; end
try epoch_end = ps.epoch_end; catch epoch_end = ''; end
try responses = ps.responses; catch responses = 0; end
try baseline = ps.baseline; catch baseline = 2; end % in seconds
try use_filtered = ps.use_filtered; catch use_filtered = 0; end
try ep_per_fig = ps.ep_per_fig; catch ep_per_fig = 3; end
try pargs = ps.pargs; catch pargs = {'-dpsc','-append'}; end
try export_stub = ps.export_stub; catch export_stub = ''; end

% check required datasets, based upon steps taken
if EXTRACT_DATA
  check_vars = {'physio_paths','pres_paths'};
  check_required_vars;
end
if CALCULATE_STATS && ~EXTRACT_DATA
  check_vars = {{'gsr_epochs','cardiac_epochs'}};
  check_required_vars;
end

% get scr analysis parameters
try scr.zscore = ps.scr.zscore; catch scr.zscore = ''; end
try scr.int_begin = ps.scr.int_begin; catch scr.int_begin = 2; end % in sec
try scr.int_end = ps.scr.int_end; catch scr.int_end = 6; end % in seconds
try scr.pk.thresh = ps.scr.find_peaks.thresh; catch scr.pk.thresh = 0; end
try scr.pk.peak_heights = ps.find_peaks.peak_heights;
  catch scr.pk.peak_heights.calcFrom = 'left';
        scr.pk.peak_heights.transformation = 'none';
end
try scr.pk.peak_height_window = ps.scr.find_peaks.peak_height_window;
  catch scr.pk.peak_height_window = 0.6; end % in seconds
try scr.pk.peak_height_thresh = ps.scr.find_peaks.peak_height_thresh;
  catch scr.pk.peak_height_thresh = 0.02; end
try keep_baseline = ps.scr.keep_baseline; catch keep_baseline = 1; end
  
% verify integral beginning and ending points
if scr.int_begin > scr.int_end
  error(['must provide valid integral beginning and ending points, in '...
      'seconds\n']);
end

% get pulse analysis parameters
try pulse.pk = ps.pulse.find_peaks;
catch pulse.pk.thresh = 0.5; pulse.pk.peak_height_window = 0.3; end

% get types of responses to include in output data
if isstruct(responses) && (isnumeric(epoch_start) && ...
        (length(epoch_start) < 1 || epoch_start)) || ...
        (isstruct(epoch_start) && isfield(epoch_start,'filter'))
  resp_names = {responses.type};
  nresp = length(resp_names);
else
  resp_names = {};
  nresp = 0;
end

% check channels against what is supported by this script
supported = {'scr','gsr','cardiac','pulse'};
nchanreq = length(chans_req);
channels = {};

% weed out unsupported channels
for ic=1:nchanreq
  if ~isempty(strmatch(chans_req{ic},supported))
    channels = [channels chans_req{ic}];
  else
    warning('channel %s not supported, not being used\n',chans_req{ic});
  end
end

% get number of supported channels that have been requested
nchan = length(channels);
if ~nchan, error('no supported channels provided'), end

eeglab('initpaths');

% init default vars used in extracted datasets
xvars = {'subject_id','session','ensemble_id','run','trial',...
    'signal','srate','peakidxs'};

if EXTRACT_DATA
  lvars = {xvars{:},resp_names{:}};
  
  if (~isempty(strmatch('gsr',channels)) || ...
          ~isempty(strmatch('scr',channels)))
    % initialize a gsr summary data output structure
    outdata.vars = [outdata.vars 'gsr_epochs'];
    gsre_idx = length(outdata.vars);
    outdata.data{gsre_idx} = ensemble_init_data_struct();
    outdata.data{gsre_idx}.type = 'gsr_epochs';
    outdata.data{gsre_idx}.vars = lvars;
    gsre_cols = set_var_col_const(outdata.data{gsre_idx}.vars);
    outdata.data{gsre_idx}.data{gsre_cols.subject_id} = {};
    outdata.data{gsre_idx}.data{gsre_cols.session} = [];
    outdata.data{gsre_idx}.data{gsre_cols.ensemble_id} = [];
    outdata.data{gsre_idx}.data{gsre_cols.run} = [];
    outdata.data{gsre_idx}.data{gsre_cols.trial} = [];
    outdata.data{gsre_idx}.data{gsre_cols.signal} = {};
    outdata.data{gsre_idx}.data{gsre_cols.srate} = [];
    outdata.data{gsre_idx}.data{gsre_cols.peakidxs} = {};
    for in=1:length(resp_names)
      outdata.data{gsre_idx}.data{gsre_cols.(resp_names{in})} = {};
    end
  end

  if ~isempty(strmatch('cardiac',channels)) || ...
          ~isempty(strmatch('pulse',channels))
    % initialize a cardiac summary data output structure
    outdata.vars = [outdata.vars 'cardiac_epochs'];
    carde_idx = length(outdata.vars);
    outdata.data{carde_idx} = ensemble_init_data_struct();
    outdata.data{carde_idx}.type = 'cardiac_epochs';
    outdata.data{carde_idx}.vars = lvars;
    carde_cols = set_var_col_const(outdata.data{carde_idx}.vars);
    outdata.data{carde_idx}.data{carde_cols.subject_id} = {};
    outdata.data{carde_idx}.data{carde_cols.session} = [];
    outdata.data{carde_idx}.data{carde_cols.ensemble_id} = [];
    outdata.data{carde_idx}.data{carde_cols.run} = [];
    outdata.data{carde_idx}.data{carde_cols.trial} = [];
    outdata.data{carde_idx}.data{carde_cols.signal} = {};
    outdata.data{carde_idx}.data{carde_cols.srate} = [];
    outdata.data{carde_idx}.data{carde_cols.peakidxs} = {};
    for in=1:length(resp_names)
      outdata.data{carde_idx}.data{carde_cols.(resp_names{in})} = {};
    end
  end
end

% % % extract data?
if EXTRACT_DATA

  %
  % START OF THE SUBJECT LOOP
  %

  for isub=1:nsub_proc
    subid = sinfo(isub).id;
    msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
    r = update_report(r,msg);

    % Determine number of sessions for this subject
    nsess = length(sinfo(isub).sessinfo);
  
    %
    % START OF THE SESSION LOOP
    %
  
    for isess = 1:nsess
      sess = sinfo(isub).sessinfo(isess);
    
      if isfield(sess,'use_session') && ~sess.use_session
        msg = sprintf('\t\t\tSkipping session %d\n', isess);
        r = update_report(r,msg);
        continue
      elseif ~isfield(sess,'physio') || ~isstruct(sess.physio)
        msg = sprintf(['\t\t\tNo Physio sinfo for session %d, '...
            'SKIPPING!\n'],isess);
        r = update_report(r,msg);
        continue
      end

      pfilt.include.all.subject_id = {subid};
      pfilt.include.all.session = isess;
      pdata = ensemble_filter(pres_paths,pfilt);
    
      if isempty(pdata.data{prescol.path}) || ...
              ~exist(pdata.data{prescol.path}{1},'file')
        error('can not find presentation data for subject %s, session %d',...
            subid,isess);
      else
        presdata = load(pdata.data{prescol.path}{1});
        presdcol = set_var_col_const(presdata.vars);
      end
    
      nruns = length(sess.use_epi_runs);
      for irun=1:nruns
      
        rnum = sess.use_epi_runs(irun);
        msg = sprintf('\t\tPROCESSING RUN %d (%d/%d)\n',rnum,irun,nruns);
        r = update_report(r,msg);
      
        lpresfilt.include.all.SUB_ID = {subid};
        lpresfilt.include.all.RUN = rnum;
      
        % get presentation data for this run
        lpres = ensemble_filter(presdata,lpresfilt);
      
        lresp = cell(nresp,1);
        for ir=1:nresp
          lrespfilt = responses(ir).event;
          lrespfilt.include.all.RUN = rnum;
          lrespdata = ensemble_filter(lpres,lrespfilt);
        
          nlresp = length(lrespdata.data{presdcol.RESP_CODE});
          lresp{ir} = cell(nlresp,1);
          for iresp=1:nlresp
            lresp{ir}{iresp} = lrespdata.data{presdcol.RESP_CODE}{iresp}{1};
          end
        end

        lphysfilt = pfilt;
        lphysfilt.include.all.run = rnum;
        if use_filtered
          lphysfilt.include.all.path = {'.*filtered.*'};
        end
      
        lphys = ensemble_filter(physio_paths,lphysfilt);
        if isempty(lphys.data{physcol.path}) || ...
                ~exist(lphys.data{physcol.path}{1},'file')
          warning('no physio file found for %s, sess %d, run %d\n',...
              subid,isess,rnum);
          continue
        end
      
        % load data
        [fp,fn,fx] = fileparts(lphys.data{physcol.path}{1});
        EEG = pop_loadset('filename',sprintf('%s%s',fn,fx),'filepath',fp);
      
        % set up epochs using presentation data
        if (isnumeric(epoch_start) && (length(epoch_start) > 1 || ...
                epoch_start)) || ...
                (isstruct(epoch_start) && isfield(epoch_start,'filter'))

          %%%% find beginnings of epochs
          % FIXME: add ability to handle multiple epoch definitions, for
          % instance a def for rest period and a def for task period
          if isnumeric(epoch_start) && length(epoch_start) > 1
            % epoch times provided in a vector, assumed in seconds
            estimes = epoch_start;
            ns = length(estimes);
          elseif isstruct(epoch_start) && isfield(epoch_start,'filter')
            % filter criteria for presentation data provided
            % events matching these criteria will be taken as epoch starts
            esf = epoch_start.filter;
            esdata = ensemble_filter(lpres,esf);
            if isempty(esdata.data{presdcol.EVENT_TYPE})
              warning(['epoch start data missing, subject %s, session '...
                  '%d, run %d\n'],subid,isess,rnum);
              continue
            end
            estimes = esdata.data{presdcol.RUN_REL_TIME}/1000;
            ns = length(estimes);
          else
            warning(['no epoch start information, skipping run %d, sub '...
                '%s, session %d\n'],rnum,subid,isess);
            continue
          end
        
          %%%% find ends of epochs
          if isnumeric(epoch_end) && length(epoch_end) == 1
            % epoch time provided, assumed in seconds
            etime = epoch_end;
          elseif isnumeric(epoch_end) && length(epoch_end) > 1
            % epoch end times provided in a vector, assumed in seconds
            % will later be compared to epoch start times
            eetimes = epoch_end;
            ne = length(eetimes);
          elseif isstruct(epoch_end) && isfield(epoch_end,'filter')
            % filter criteria for presentation data provided
            % events matching these criteria will be taken as epoch ends
            eef = epoch_end.filter;
            eedata = ensemble_filter(lpres,eef);
            if isempty(eedata.data{presdcol.EVENT_TYPE})
              warning(['epoch end data missing, subject %s, session '...
                  '%d, run %d\n'],subid,isess,rnum);
              continue
            end
            eetimes = eedata.data{presdcol.RUN_REL_TIME}/1000;
            ne = length(eetimes);
          else
            warning(['no epoch end information provided for sub %s, sess '...
                '%d, run %d, calculating as minimum distance between '...
                'epoch onsets\n'],subid,isess,rnum);
            etime = min(diff(estimes));
          end
          
          % account for raw signal baseline period
          estimes = estimes - EEG.xmin;
          eetimes = eetimes - EEG.xmin;

          if ~exist('etime','var') && ne ~= ns || any (estimes > eetimes)
            % end time not specified, 
            warning(['outlier starting or ending events, or events not '...
                'aligned, sub %s, session %d, run %d\n'],subid,isess,rnum);
            continue
            % FIXME: should have some logic here to at least attempt to
            % reconcile estimes and eetimes ... rules such as "no end time
            % before the first start, no start time before the first end", as
            % well as "line up each start/end, 1 to 1, pick out outliers, and
            % reconcile pairs whose end is before it's given start          
          end

          % find epoch length, if not explicitly provided
          if ~exist('etime','var') 
            % make sure to remove epochs that exceed the samples in EEG
            ee_end_idxs = eetimes > (EEG.pnts/EEG.srate);
            if any(ee_end_idxs)
              warning('%d epochs out of bounds, removing\n',length(ee_end_idxs));
            end
            eetimes(ee_end_idxs) = []; % remove exceeding epochs from eetimes
            estimes(ee_end_idxs) = []; % remove exceeding epochs from estimes
            dtimes = eetimes - estimes;
            if isstruct(epoch_end) && isfield(epoch_end,'minimum_end') ...
                    && epoch_end.minimum_end
              etime = min(dtimes);
            else
              etime = max(dtimes);
            end
          end

          % double check that all epochs end before the end of the data set
          badtimes = (estimes+etime) > (EEG.pnts/EEG.srate);
          if any(badtimes)
            ee_end_idxs = find(badtimes);
            warning('%d epochs out of bounds, removing\n',length(ee_end_idxs));
            eetimes(ee_end_idxs) = []; % remove exceeding epochs from eetimes
            estimes(ee_end_idxs) = []; % remove exceeding epochs from estimes
          end
          ns = length(estimes);
          ne = length(eetimes);

          % import event onsets
          EEG = pop_importevent(EEG,'append','no','event',[ones(ns,1) ...
              estimes],'fields',{'type','latency'},'timeunit',1);

          % epoch the data based upon event onsets and the largest event duration 
          EEG = pop_epoch(EEG,[],[-baseline etime],'eventindices',1:ns);

          if EEG.trials ~= ns
            error('one or more epochs was lost, sub %s, sess %d, run %d\n',...
                subid,isess,rnum);
          end
          
        else
          ns = 1;
        end % if (isnumeric(epoch_start
        
        % set sampling rate
        scr.pk.Fs = EEG.srate;
        pulse.pk.Fs = EEG.srate;

        % set up a figure file for tsub/sess/run
        if GRAPH_EPOCHS && SAVE_EPOCHS
          figfname = fullfile(fp,sprintf('%s_epochs.ps',fn));
        end
        
        for ic=1:nchan
          c = channels{ic};

          baseline_end = (baseline*EEG.srate)+1; % next sample after end of baseline

          % extract channel location
          cidx = strmatch(c,{EEG.chanlocs.labels});
          if isempty(cidx)
            warning(['no channel label found for %s, skipping for '...
                'sub %s, session %d, run %d'],c,subid,isess,rnum);
            continue
          end

          switch c
            case {'gsr','scr'}
            
              % extract all epochs for this channel
              epochs = squeeze(EEG.data(cidx,:,:));
              
              sepoch = size(epochs);
              if (sepoch(1) == ns && sepoch(2) ~= ns)
                epochs = epochs';
                sepoch = size(epochs);
              end

              % subtract baseline for each epoch
              for iep=1:ns
                mbase = mean(epochs(1:(baseline*EEG.srate),iep));
                pksig = epochs(:,iep) - mbase;
                
                if ~keep_baseline && baseline
                  pksig(1:baseline*EEG.srate) = [];
                end
                
                % get peaks
                [pidxs,heights,bidxs] = find_peaks(pksig,scr.pk);

                if keep_baseline && baseline
                  % remove peaks whose troughs occur before 1 sec after
                  % stimulus onset
                  rmbase = bidxs < (baseline+1)*EEG.srate;
                  pidxs(rmbase) = [];
                  heights(rmbase) = [];
                  bidxs(rmbase) = [];
                end
                
                % graph??
                if GRAPH_EPOCHS
                  m = mod(iep-1,ep_per_fig);
                  if m == 0, figure(); end
                  subplot(ep_per_fig,1,m+1);
                  plot(pksig);
                  title(sprintf(['%s, sess %d, run %d, signal %s, epoch'...
                      ' %d'],subid,isess,rnum,c,iep));
                  hold on;
                  for ipk = 1:length(pidxs)
                    plot(pidxs(ipk),pksig(pidxs(ipk)),'g*');
                  end
                  hold off;
                  
                  if SAVE_EPOCHS && (m+1) == ep_per_fig
                    % save graphed epochs to a file
                    print(pargs{:},figfname);
                  end
                end
                
                % save signal out to gsr_epochs
                outdata.data{gsre_idx} = ensemble_add_data_struct_row(...
                    outdata.data{gsre_idx},'subject_id',subid,'session',...
                    isess,'ensemble_id',sess.ensemble_id,'run',rnum,...
                    'trial',iep,'signal',{pksig},'srate',EEG.srate,...
                    'peakidxs',pidxs);
              end % for iep=1:ns
              
              % save final figure to file if not already saved
              if SAVE_EPOCHS && (m+1) ~= ep_per_fig
                print(pargs{:},figfname);
              end

            case {'cardiac','pulse'}
                
              % iterate over epochs
              for iep=1:ns
                % find peaks
                pksig = EEG.data(cidx,baseline_end:end,iep);
                pidxs = find_peaks(pksig,pulse.pk);

                % save signal out to cardiac_epochs
                outdata.data{carde_idx} = ensemble_add_data_struct_row(...
                    outdata.data{carde_idx},'subject_id',subid,'session',...
                    isess,'ensemble_id',sess.ensemble_id,'run',rnum,...
                    'trial',iep,'signal',{pksig'},'srate',EEG.srate,...
                    'peakidxs',pidxs);

                % graph??
                if GRAPH_EPOCHS
                  m = mod(iep-1,ep_per_fig);
                  if m == 0, figure(); end
                  subplot(ep_per_fig,1,m+1);
                  plot(pksig);
                  title(sprintf(['%s, sess %d, run %d, signal %s, epoch'...
                      ' %d'],subid,isess,rnum,c,iep));
                  hold on;
                  for ipk = 1:length(pidxs)
                    plot(pidxs(ipk),pksig(pidxs(ipk)),'g*');
                  end
                  hold off;

                  if SAVE_EPOCHS && (m+1) == ep_per_fig
                    % save graphed epochs to a file
                    print(pargs{:},figfname);
                  end
                end
              end % for iep=1:ns
              
              % save final figure to file if not already saved
              if SAVE_EPOCHS && (m+1) ~= ep_per_fig
                print(pargs{:},figfname);
              end

            otherwise
              warning('unknown channel %s',c);
              
          end % switch c

          % add response information
          for iresp = 1:length(resp_names)
            % get all responses
            rname = resp_names{iresp};
            respdata = lresp{iresp};
            nrespd = length(respdata);
            % check against length of events (ns)
            if nrespd < ns
              warning(['data for response %s of different length (%d) '...
                  'than number of epochs (%d), these are not being '...
                  'added to the output\n'],rname,nrespd,ns);
              continue
            elseif ns < nrespd
              warning(['more responses (%d) than epochs (%d), assume '...
                  'that epochs were dropped from the end, and dropping '...
                  'extra responses'],nrespd,ns);
              respdata(ns+1:end) = [];
            end

            % get indices for this channel and this response
            switch c
              case {'gsr','scr'}
                oidx = gsre_idx;
                ridx = gsre_cols.(rname);
              case {'cardiac','pulse'}
                oidx = carde_idx;
                ridx = carde_cols.(rname);
              otherwise
                warning('unknown channel type (%s)\n',c);
                continue
            end

            % save into outdata
            outdata.data{oidx}.data{ridx}(end-ns+1:end) = respdata;
            
          end % for iresp=1:length(resp_names
        end % for ic=1:nchan
      end % for irun=1:
    end % for isess = 1:

    if (GRAPH_EPOCHS && EPOCH_BY_RESP)
      pfilt = struct();
      pfilt.include.all.subject_id = {subid};
      pfilt.include.all.path_type = {'physio_outdir'};
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        physdir = lpathdata.data{pcol.path}{1};
        check_dir(physdir);
      else
        warning(['could not find physio_outdir in pathdata, not graphing '...
            'epochs by response for subject %s'],subid);
        continue
      end

      if GRAPH_EPOCHS && EPOCH_BY_RESP
        figfname = fullfile(physdir,sprintf('%s_epochs_by_resp.ps',subid));
        % graph epochs by response type
        for iresp=1:length(resp_names)
          rname = resp_names{iresp};
          for ic=1:nchan
            c = channels{ic};
            % get channel and response indexes
            switch c
              case {'gsr','scr'}
                eidx = gsre_idx;
                ridx = gsre_cols.(rname);
                sidx = gsre_cols.signal;
              case {'cardiac','pulse'}
                eidx = carde_idx;
                ridx = carde_cols.(rname);
                sidx = carde_cols.signal;
              otherwise
                warning('unknown channel type (%s)\n',c);
                continue
            end % switch c

            % get mask matrix
            if iscell(outdata.data{eidx}.data{ridx}) && ...
                    all(cellfun(@isempty,outdata.data{eidx}.data{ridx}))
              warning('no responses for response_name %s\n',rname);
              continue
            end

            [rm,ur] = make_mask_mtx(outdata.data{eidx}.data{ridx});
        
            for ir=1:length(ur)
              figure();
              allsig = [];
              hold on;
              lidxs = find(rm(:,ir));
              for i=1:length(lidxs)
                li=lidxs(i);
                sig=outdata.data{eidx}.data{sidx}{li};              
                if isempty(sig), continue, end
                if ~isempty(strmatch(c,{'cardiac','pulse'}))
                  pidxs=find_peaks(sig,pulse.pk);
                  if isempty(pidxs), continue, end
                  sig=60./(diff(pidxs)/EEG.srate);
                end
                if size(sig,1) > size(sig,2), sig = sig'; end
                plot(sig,'g');
                ssdiff = size(allsig,2) - size(sig,2);
                if ssdiff > 0
                  % allsig is longer, extend sig with nans at the end
                  sig = [sig nan(1,ssdiff)];
                elseif ssdiff < 0
                  % sig is longer, extend allsig with nans at the end
                  allsig = [allsig nan(size(allsig,1),-ssdiff)];
                end
                allsig = [allsig; sig];
              end % for i=1:length(lidxs
              plot(nanmean(allsig),'k');
              hold off;
              lr = ur{ir};
              if isnumeric(lresp), lresp = num2str(lresp); end
              title(sprintf('signal: %s, response: %s, %s',c,r,lr));
          
              if SAVE_EPOCHS, print(pargs{:},figfname); end
            end % for ir=1:length(ur
          end % for iresp=1:length(resp_names
        end % for ic=1:nchan
      end % if GRAPH_EPOCHS && EPOCH_BY_RESP
    end % if (GRAPH_EPOCHS && EPOCH_BY_RESP) ||

  end % for isub=1:

end % if EXTRACT_DATA

if CALCULATE_STATS

  % get gsr data?
  if exist('gsr_epochs','var')
    gdata = gsr_epochs;
  elseif exist('gsre_idx','var')
    gdata = outdata.data{gsre_idx};
  else
    warning('no gsr epoched data found');
  end
  
  % get cardiac data?
  if exist('cardiac_epochs','var')
    cdata = cardiac_epochs;
  elseif exist('carde_idx','var')
    cdata = outdata.data{carde_idx};
  else
    warning('no cardiac epoched data found');
  end
  
  if exist('gdata','var')
      
    % add response information
    ridxs = ~ismember(gdata.vars,xvars);
    rvars = {gdata.vars{ridxs}};
    if ~iscell(rvars), rvars = {rvars}; end

    % initialize a gsr summary data output structure
    outdata.vars = [outdata.vars 'gsr_summ'];
    gsr_idx = length(outdata.vars);
    outdata.data{gsr_idx} = ensemble_init_data_struct();
    outdata.data{gsr_idx}.type = 'gsr_summ';
    outdata.data{gsr_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','integral','scl','ttl_integral','ttl_scl',...
        'npeaks','mean_peak_height','peak_height_variance',...
        'first_peak_time','first_peak_height',rvars{:}};
    gsr_cols = set_var_col_const(outdata.data{gsr_idx}.vars);
    outdata.data{gsr_idx}.data{gsr_cols.subject_id} = {};
    outdata.data{gsr_idx}.data{gsr_cols.session} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ensemble_id} = [];
    outdata.data{gsr_idx}.data{gsr_cols.run} = [];
    outdata.data{gsr_idx}.data{gsr_cols.trial} = [];
    outdata.data{gsr_idx}.data{gsr_cols.integral} = [];
    outdata.data{gsr_idx}.data{gsr_cols.scl} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ttl_integral} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ttl_scl} = [];
    outdata.data{gsr_idx}.data{gsr_cols.npeaks} = [];
    outdata.data{gsr_idx}.data{gsr_cols.mean_peak_height} = [];
    outdata.data{gsr_idx}.data{gsr_cols.peak_height_variance} = [];
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_time} = [];
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_height} = [];
    for in=1:length(rvars)
      outdata.data{gsr_idx}.data{gsr_cols.(rvars{in})} = {};
    end

    % initialize gsr vars
    gcol = set_var_col_const(gdata.vars);
    ng = length(gdata.data{1});
    
    % calculate stats for each epoch
    for ie=1:ng
      sig   = gdata.data{gcol.signal}{ie};
      srate = gdata.data{gcol.srate}(ie);
      pidxs = gdata.data{gcol.peakidxs}{ie};
      if iscell(pidxs), pidxs = pidxs{1}; end
      if ~isempty(epoch_end) && isnumeric(epoch_end)
        end_samp = epoch_end*srate;
        sig = sig(1:end_samp);
        pidxs(pidxs > end_samp) = [];
      end
      if isempty(pidxs)
        pidxs = [];
        hghts = [];
        mpkhght = NaN;
        pkhghtv = NaN;
        fpktime = NaN;
        fpkhght = NaN;
      else
        hghts = peak_heights('signal',sig,'params',scr.pk.peak_heights,...
            'peakIdxs',pidxs);
        mpkhght = mean(hghts);
        pkhghtv = var(hghts);
        fpktime = pidxs(1)/srate;
        fpkhght = hghts(1);
      end
      npeaks = length(pidxs);
      
      baseline_end = baseline*srate;

      % calculate integral within integral window defined by
      % int_begin and int_end
      startint = baseline_end + scr.int_begin*srate + 1;
      endint   = min([(baseline_end + scr.int_end*srate) length(sig)]);

      integral = trapz(sig(startint:endint));
      mscl = mean(sig(startint:endint));
      ttl_int = trapz(sig(baseline_end+1:end));
      ttlmscl = mean(sig(baseline_end+1:end));

      % save to outdata
      outdata.data{gsr_idx} = ensemble_add_data_struct_row(...
          outdata.data{gsr_idx},...
          'subject_id',gdata.data{gcol.subject_id}{ie},...
          'session',gdata.data{gcol.session}(ie),...
          'ensemble_id',gdata.data{gcol.ensemble_id}(ie),...
          'run',gdata.data{gcol.run}(ie),...
          'trial',gdata.data{gcol.trial}(ie),...
          'integral',integral(:),'scl',mscl(:),'ttl_integral',...
          ttl_int(:),'ttl_scl',ttlmscl(:),...
          'npeaks',npeaks,'mean_peak_height',mpkhght,...
          'peak_height_variance',pkhghtv,...
          'first_peak_time',fpktime,'first_peak_height',fpkhght);

    end % for ie=1:ng
    
    % pass on response information
    for ir = 1:length(rvars)
      % get column indices
      ri = gcol.(rvars{ir});
      ro = gsr_cols.(rvars{ir});

      % save into outdata
      outdata.data{gsr_idx}.data{ro} = gdata.data{ri};
    end

  end %  if exist('gdata
  
  if exist('cdata','var')
      
    % add response information
    ridxs = ~ismember(cdata.vars,xvars);
    rvars = {cdata.vars{ridxs}};
    if ~iscell(rvars), rvars = {rvars}; end

    % initialize a cardiac summary data output structure
    outdata.vars = [outdata.vars 'cardiac_summ'];
    card_idx = length(outdata.vars);
    outdata.data{card_idx} = ensemble_init_data_struct();
    outdata.data{card_idx}.type = 'cardiac_summ';
    outdata.data{card_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','heart_rate','hr_variance','hr_slope',rvars{:}};
    card_cols = set_var_col_const(outdata.data{card_idx}.vars);
    outdata.data{card_idx}.data{card_cols.subject_id} = {};
    outdata.data{card_idx}.data{card_cols.session} = [];
    outdata.data{card_idx}.data{card_cols.ensemble_id} = [];
    outdata.data{card_idx}.data{card_cols.run} = [];
    outdata.data{card_idx}.data{card_cols.trial} = [];
    outdata.data{card_idx}.data{card_cols.heart_rate} = [];
    outdata.data{card_idx}.data{card_cols.hr_variance} = [];
    outdata.data{card_idx}.data{card_cols.hr_slope} = [];
    for in=1:length(rvars)
      outdata.data{card_idx}.data{card_cols.(rvars{in})} = {};
    end

    % initialize cardiac vars
    ccol = set_var_col_const(cdata.vars);
    nc = length(cdata.data{1});
    
    % iterate over epochs
    for iep=1:nc
      srate = cdata.data{ccol.srate}(iep);
      pidxs = cdata.data{ccol.peakidxs}{iep};
      if iscell(pidxs), pidxs = pidxs{1}; end
      if ~isempty(epoch_end) && isnumeric(epoch_end)
        end_samp = epoch_end*srate;
        sig = sig(1:end_samp);
        pidxs(pidxs > end_samp) = [];
      end
      if isempty(pidxs), pidxs = nan; end

      dtimes = diff(pidxs);

      % remove extremely small dtimes
      dtimes(dtimes < mean(dtimes)-3*std(dtimes)) = [];

      % heart rate (in beats per minute) and heart rate variance
      hr_bpm = 1/(mean(dtimes)/srate)*60;
      hr_var = var(dtimes);

      % slope of change in heart rate
      zlength = zscore(1:length(dtimes))';
      ztimes = zscore(dtimes);
      hr_slope = regress(ztimes,zlength);

      % save to outdata
      outdata.data{card_idx} = ensemble_add_data_struct_row(...
          outdata.data{card_idx},...
          'subject_id', cdata.data{ccol.subject_id}{iep},...
          'session',    cdata.data{ccol.session}(iep),...
          'ensemble_id',cdata.data{ccol.ensemble_id}(iep),...
          'run',        cdata.data{ccol.run}(iep),...
          'trial',      cdata.data{ccol.trial}(iep),...
          'heart_rate',hr_bpm,'hr_variance',hr_var,'hr_slope',hr_slope);

    end % for iep=1:nc

    
    % pass on response information
    for ir = 1:length(rvars)
      % get column indices
      ri = ccol.(rvars{ir});
      ro = card_cols.(rvars{ir});

      % save into outdata
      outdata.data{card_idx}.data{ro} = cdata.data{ri};
    end
    
  end % if exist('cdata

  if SAVE_DATA
    for ic=1:nchan
      c = channels{ic};
      % get channel idxs
      switch c
        case {'gsr','scr'}
          oidx = gsr_idx;
          ocol = gsr_cols;
        case {'cardiac','pulse'}
          oidx = card_idx;
          ocol = card_cols;
        otherwise
          warning('unknown channel type (%s)\n',c);
      end % switch c

      % find unique subjects in calculated data
      us = unique(outdata.data{oidx}.data{ocol.subject_id});
      ns = length(us);
      
      % iterate over subjects
      for is=1:ns
        subid = us{is};
        
        % get output data for this subject
        pfilt = struct();
        pfilt.include.all.subject_id = {subid};
        odata = ensemble_filter(outdata.data{oidx},pfilt);
        
        % get path data for this subject
        pfilt.include.all.path_type = {'physio_outdir'};
        lpathdata = ensemble_filter(pathdata,pfilt);
        if isempty(lpathdata.data{1})
          warning(['could not find physio_outdir in pathdata, not '...
              'saving data for subject %s'],subid);
        elseif isempty(odata.data{1})
          warning('could not find output data for subject %s, SKIPPING',...
              subid);
        else
          physdir = lpathdata.data{pcol.path}{1};
          check_dir(physdir);

          xsp.export = export;
          xsp.export.fname = fullfile(physdir,...
              sprintf('%s_%s_%sexport.txt',subid,odata.type,export_stub));
          xsp.sas = sas;
          xsp.sas.fname = fullfile(physdir,sprintf('%s_%s_%sexport.sas',...
              subid,odata.type,export_stub));
          xsp.sas.libname=sprintf('s%s_%s%s',subid,odata.type,export_stub);
          ensemble_export_sastxt(odata,xsp);
        end % if ~isempty(lpathdata.data{1
      end % for is=1:ns
    end % for ic=1:nchan
  end % if SAVE_DATA

end % if CALCULATE_STATS
