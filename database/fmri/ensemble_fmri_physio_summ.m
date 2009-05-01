function outdata = ensemble_fmri_physio_summ(indata,defs)

% calculates summary statistics for physio data collected during fMRI
% 
% LOADING PHYSIO DATA FILES
%   defs.physio_summ.use_filtered (default 1)
% 
% EPOCHING
%   defs.physio_summ.epoch_start (default 0) - if a struct, with sub-field
%       'filter', this is treated as a filter criteria for the included
%       presentation data, to filter for events to be treated as epoch
%       onsets. If a numeric vector, treated as a vector of epochs onsets
%       in seconds.
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
%       will be used for all epochs.
%   defs.physio_summ.lead (default 2) - pre-epoch baseline to be added to
%       epoch when calculating and extracting epochs. this value is assumed
%       to be in 'seconds'
%   measures
%   exporting
% 
% SUMMARY DATA CALCULATION
%   defs.physio_summ.channels - channel names, as they appear in your
%       eeglab chanlocs labels, that you wish to calculate summary
%       information for.
%   defs.physio_summ.scr - settings for SCR/GSR data summary
%       .spc_integral.begin (default 2) - time, in seconds, to start
%           calculation of the specific SCR integral
%       .spc_integral.end (default 6) - time, in seconds, to end
%           calculation of the specific SCR integral
%   	.proxlim (default 0.6) - limit of proximity of peak heights, in
%           seconds. for each peak found, if there are other peaks within
%           +/- proxlim seconds of that peak, all but the highest peaks
%           will be discarded.
%   	.find_peaks.thresh (default 0) - proportion of total signal
%           range below which peaks will not be considered. for example,
%           for .thresh = 0.5, peaks will only be returned if their height 
%           is more than half of the total signal range
%   	.peak_heights.calcFrom (default 'left') - parameter passed to
%           'peak_heights', indicating from which side of a peak to find the
%           trough used to calculate peak height
%       .peak_heights.transformation (default 'none') - parameter passed to
%           'peak_heights' to indicate the transformation, if any, to
%           perform on peak height values before returning them.
%   	.height_thresh (default 0.02) - peak heights below this value will
%           be rejected
%       .zscore (default '') - if set to 1, each scr epoch will be zscored
%           within-epoch, before being analyzed
%   defs.physio_summ.pulse - settings for pulse ox. data summary
%   	.find_peaks.thresh (default 0.5) - proportion of total signal
%           range below which peaks will not be considered. for example,
%           for .thresh = 0.5, peaks will only be returned if their height 
%           is more than half of the total signal range
% 
% RESPONSE MATCHING FOR EPOCHS
%   presentation data must be included, in the form of a pres_paths data
%       struct (created by presentation_parse) containing paths to
%       presentation data that has been parsed and saved to disk for the
%       given participant/session.
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
% DATA OUTPUT
% 
%   defs.physio_summ.GRAPH_EPOCHS (default 0)
%   defs.physio_summ.EPOCHS_BY_RESP (default 0)
%   defs.physio_summ.SAVE_EPOCHS (default 1)
%   defs.physio_summ.SAVE_DATA (default 0)
%   defs.physio_summ.export (default 0)
%   defs.physio_summ.sas (default 0)
% 
% FIXME - should be re-written, to include a mechanism by which channels to
% analyze are identified by EEG.chanlocs name as sub-fields of the params
% struct, and under that, sub-fields of each channel should be named
% exactly as the names of functions used to either filter or report the
% data. Then, sub-fields of those fields should be params that get sent to
% those functions as the second argument, the first argument being the
% data output of the previous function (in may or all cases, an eeglab
% struct?
% 
% FB 2009.04.21

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'fmri_physio_summ';  % Identify the type of this reporting instance
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
    end
  end
end

% check for required vars
check_vars = {'sinfo','pathdata','physio_paths','pres_paths'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && ~isempty(pathdata.data{1})
    if length(nsub_proc) == 1
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
try export = ps.export; catch export.write2file = 0; end
try sas = ps.sas; catch sas.libsave = 0; end
try epoch_start = ps.epoch_start; catch epoch_start = ''; end
try epoch_end = ps.epoch_end; catch epoch_end = ''; end
try responses = ps.responses; catch responses = 0; end
try lead = ps.lead; catch lead = 2; end % in seconds
try use_filtered = ps.use_filtered; catch use_filtered = 0; end
try ep_per_fig = ps.ep_per_fig; catch ep_per_fig = 3; end
try pargs = ps.pargs; catch pargs = {'-dpsc','-append'}; end

% get scr analysis parameters
try scr.zscore = ps.scr.zscore; catch scr.zscore = ''; end
try scr.int_begin = ps.scr.int_begin; catch scr.int_begin = 2; end % in sec
try scr.int_end = ps.scr.int_end; catch scr.int_end = 6; end % in seconds
try scr.proxlim = ps.scr.proxlim; catch scr.proxlim = 0.6; end % in seconds
try scr.pkh.calcFrom = ps.peak_heights.calcFrom;
  catch scr.pkh.calcFrom = 'left'; end
try scr.pkh.transformation = ps.peak_heights.transformation;
  catch scr.pkh.transformation = 'none'; end
try scr.pk.thresh = ps.scr.find_peaks.thresh; catch scr.pk.thresh = 0; end
try scr.height_thresh = ps.scr.height_thresh; catch scr.height_thresh = 0.02; end

% verify integral beginning and ending points
if scr.int_begin > scr.int_end
  error(['must provide valid integral beginning and ending points, in '...
      'seconds\n']);
end

% get pulse analysis parameters
try pulse.pk.thresh = ps.pulse.find_peaks.thresh;
catch pulse.pk.thresh = 0.5; end

% get types of responses to include in output data
if isstruct(responses)
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

if ~isempty(strmatch('gsr',channels)) || ~isempty(strmatch('scr',channels))
    % initialize a gsr summary data output structure
    outdata.vars = [outdata.vars 'gsr_summ'];
    gsr_idx = length(outdata.vars);
    outdata.data{gsr_idx} = ensemble_init_data_struct();
    outdata.data{gsr_idx}.type = 'gsr_summ';
    outdata.data{gsr_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','integral','mean_activation','ttl_integral',...
        'ttl_mean_activation','npeaks','mean_peak_height',...
        'peak_height_variance','first_peak_time','first_peak_height',...
        resp_names{:}};
    gsr_cols = set_var_col_const(outdata.data{gsr_idx}.vars);
    outdata.data{gsr_idx}.data{gsr_cols.subject_id} = {};
    outdata.data{gsr_idx}.data{gsr_cols.session} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ensemble_id} = [];
    outdata.data{gsr_idx}.data{gsr_cols.run} = [];
    outdata.data{gsr_idx}.data{gsr_cols.trial} = [];
    outdata.data{gsr_idx}.data{gsr_cols.integral} = [];
    outdata.data{gsr_idx}.data{gsr_cols.mean_activation} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ttl_integral} = [];
    outdata.data{gsr_idx}.data{gsr_cols.ttl_mean_activation} = [];
    outdata.data{gsr_idx}.data{gsr_cols.npeaks} = [];
    outdata.data{gsr_idx}.data{gsr_cols.mean_peak_height} = [];
    outdata.data{gsr_idx}.data{gsr_cols.peak_height_variance} = [];
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_time} = [];
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_height} = [];
    for in=1:length(resp_names)
      outdata.data{gsr_idx}.data{gsr_cols.(resp_names{in})} = {};
    end
    
    % initialize a gsr summary data output structure
    outdata.vars = [outdata.vars 'gsr_epochs'];
    gsre_idx = length(outdata.vars);
    outdata.data{gsre_idx} = ensemble_init_data_struct();
    outdata.data{gsre_idx}.type = 'gsr_epochs';
    outdata.data{gsre_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','signal'};
    gsre_cols = set_var_col_const(outdata.data{gsre_idx}.vars);
    outdata.data{gsre_idx}.data{gsre_cols.subject_id} = {};
    outdata.data{gsre_idx}.data{gsre_cols.session} = [];
    outdata.data{gsre_idx}.data{gsre_cols.ensemble_id} = [];
    outdata.data{gsre_idx}.data{gsre_cols.run} = [];
    outdata.data{gsre_idx}.data{gsre_cols.trial} = [];
    outdata.data{gsre_idx}.data{gsre_cols.signal} = {};
end

if ~isempty(strmatch('cardiac',channels)) || ...
        ~isempty(strmatch('pulse',channels))
    % initialize a cardiac summary data output structure
    outdata.vars = [outdata.vars 'cardiac_summ'];
    card_idx = length(outdata.vars);
    outdata.data{card_idx} = ensemble_init_data_struct();
    outdata.data{card_idx}.type = 'cardiac_summ';
    outdata.data{card_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','heart_rate','hr_variance',resp_names{:}};
    card_cols = set_var_col_const(outdata.data{card_idx}.vars);
    outdata.data{card_idx}.data{card_cols.subject_id} = {};
    outdata.data{card_idx}.data{card_cols.session} = [];
    outdata.data{card_idx}.data{card_cols.ensemble_id} = [];
    outdata.data{card_idx}.data{card_cols.run} = [];
    outdata.data{card_idx}.data{card_cols.trial} = [];
    outdata.data{card_idx}.data{card_cols.heart_rate} = [];
    outdata.data{card_idx}.data{card_cols.hr_variance} = [];
    for in=1:length(resp_names)
      outdata.data{card_idx}.data{card_cols.(resp_names{in})} = {};
    end

    % initialize a cardiac summary data output structure
    outdata.vars = [outdata.vars 'cardiac_epochs'];
    carde_idx = length(outdata.vars);
    outdata.data{carde_idx} = ensemble_init_data_struct();
    outdata.data{carde_idx}.type = 'cardiac_epochs';
    outdata.data{carde_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','signal'};
    carde_cols = set_var_col_const(outdata.data{carde_idx}.vars);
    outdata.data{carde_idx}.data{carde_cols.subject_id} = {};
    outdata.data{carde_idx}.data{carde_cols.session} = [];
    outdata.data{carde_idx}.data{carde_cols.ensemble_id} = [];
    outdata.data{carde_idx}.data{carde_cols.run} = [];
    outdata.data{carde_idx}.data{carde_cols.trial} = [];
    outdata.data{carde_idx}.data{carde_cols.signal} = {};
end

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
      if ~isempty(epoch_start) && ~isempty(epoch_end)

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
          ee_end_idxs = eetimes > EEG.pnts/EEG.srate;
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

        % import event onsets
        EEG = pop_importevent(EEG,'append','no','event',[ones(ns,1) ...
            estimes],'fields',{'type','latency'},'timeunit',1);

        % epoch the data based upon event onsets and the largest event duration 
        EEG = pop_epoch(EEG,[],[-lead etime],'eventindices',1:ns);

        if EEG.trials ~= ns
          error('one or more epochs was lost, sub %s, sess %d, run %d\n',...
              subid,isess,rnum);
        end
                
        % fill out an ns-length cell and array for session/run info
        lsubids = cell(ns,1);
        for iep=1:ns, lsubids{iep} = subid; end
        lsessid = repmat(isess,ns,1);
        lrun = repmat(rnum,ns,1);
        lesess = repmat(sess.ensemble_id,ns,1);
        
        % set up a figure file for tsub/sess/run
        if GRAPH_EPOCHS && SAVE_EPOCHS
          figfname = fullfile(fp,sprintf('%s_epochs.ps',fn));
        end
        
        for ic=1:nchan
          c = channels{ic};

          lead_end = (lead*EEG.srate)+1;

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
              
              if ~isempty(scr.zscore)
                epochs = (epochs-mean(epochs(:)))/std(epochs(:));
              end

              % subtract baseline for each epoch
              for iep=1:ns
                mbase = mean(epochs(1:(lead*EEG.srate),iep));
                epochs(:,iep) = epochs(:,iep) - mbase;
              end
          
              % calculate integral within integral window defined by
              % int_begin and int_end
              startint = lead_end + scr.int_begin*EEG.srate + 1;
              endint   = lead_end + scr.int_end*EEG.srate + 1;
              integral = trapz(squeeze(epochs(startint:endint,:)));
              integral = integral(:);
              
              % calculate mean activation across integral window
              mact = mean(squeeze(epochs(startint:endint,:)));
              mact = mact(:);
              
              % calculate integral across entire waveform
              ttl_int = trapz(squeeze(epochs(lead_end+1:end,:)));
              ttl_int = ttl_int(:);
              
              % calculate mean activation across entire waveform
              ttlmact = mean(squeeze(epochs(lead_end+1:end,:)));
              ttlmact = ttlmact(:);
              
              % init peak stats vectors
              npeaks = zeros(ns,1); % total number of peaks
              fptime = zeros(ns,1); % time of the first peak, in seconds
              fphght = zeros(ns,1); % height of the first peak, relative to last trough
              mpkhgt = zeros(ns,1); % mean peak height
              vpkhgt = zeros(ns,1); % peak height variance
              
              % iterate over epochs
              for iep=1:ns
                pksig = epochs(:,iep);
                
                % get peaks
                pidxs = find_peaks(pksig,scr.pk);

                % get peak heights
                [heights,bidxs]= peak_heights...
                    ('signal',pksig,'peakIdxs',pidxs,'params',scr.pkh);
                
                % remove peaks that do not meet the threshold criterion
                peakmask = heights > scr.height_thresh;
                pidxs(~peakmask)=[];
                heights(~peakmask)=[];
                bidxs(~peakmask)=[];

                % remove peaks whose troughs occur before 1 sec after
                % stimulus onset
                pidxs(bidxs < (lead+1)*EEG.srate) = [];
                
                % remove peaks within +/- 'proxlim' seconds of each peak
                pidxs_mod = pidxs;
                for ip=1:length(pidxs)
                    
                  % get peak window
                  peak_start=pidxs(ip)-round(EEG.srate*scr.proxlim);
                  peak_end  =pidxs(ip)+round(EEG.srate*scr.proxlim);
                  if peak_start < 0, peak_start = 0; end
                  if peak_end > EEG.pnts, peak_end = EEG.pnts; end
                  peakwindow= (peak_start:peak_end);

                  % get peaks within window
                  windowidxs=intersect(peakwindow, pidxs);
                  
                  % get peak heights, remove max peak from windowidxs
                  [windowheights]= peak_heights('signal',pksig,...
                      'peakIdxs',windowidxs,'params',scr.pkh);
                  windowidxs(windowheights==max(windowheights)) = [];

                  % remove windowidxs from pidxs_mod
                  if ~isempty(windowidxs) && ~isempty(pidxs_mod)
                    pidxs_mod(ismember(pidxs_mod,windowidxs)) = [];
                  end
                end
                pidxs = pidxs_mod;
                
                if ~isempty(pidxs)
                  % get peak heights
                  heights = peak_heights('signal',pksig,'params',...
                      scr.pkh,'peakIdxs',pidxs);

                  % save to peak stat vectors
                  npeaks(iep) = length(pidxs);
                  fptime(iep) = pidxs(1)/EEG.srate;
                  fphght(iep) = heights(1);
                  mpkhgt(iep) = mean(heights);
                  vpkhgt(iep) = var(heights);
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
                  for ipk = 1:npeaks(iep)
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
                  'trial',iep,'signal',{pksig});
              end % for iep=1:ns
              
              % save final figure to file if not already saved
              if SAVE_EPOCHS && (m+1) ~= ep_per_fig
                print(pargs{:},figfname);
              end

              % save to outdata
              outdata.data{gsr_idx} = ensemble_add_data_struct_row(...
                  outdata.data{gsr_idx},'subject_id',lsubids,'session',...
                  lsessid,'ensemble_id',lesess,'run',lrun,'trial',...
                  [1:ns]','integral',integral,'mean_activation',...
                  mact,'ttl_integral',ttl_int,'ttl_mean_activation',...
                  ttlmact,'npeaks',npeaks,'mean_peak_height',mpkhgt,...
                  'peak_height_variance',vpkhgt,'first_peak_time',...
                  fptime,'first_peak_height',fphght);

            case {'cardiac','pulse'}
                
              hr_bpm = zeros(ns,1); % avg epoch heart rate, beats per min.
              hr_var = zeros(ns,1); % heart rate variance during epoch
              
              % iterate over epochs
              for iep=1:ns
                % find peaks
                pksig = EEG.data(cidx,lead_end+1:end,iep);
                pidxs = find_peaks(pksig,pulse.pk);
                dtimes = diff(pidxs);
                
                % remove extremely small dtimes
                dtimes(dtimes < mean(dtimes)-3*std(dtimes)) = [];

                % heart rate (in beats per minute) and heart rate variance
                hr_bpm(iep) = 1/(mean(dtimes)/EEG.srate)*60;
                hr_var(iep) = var(dtimes);

                % save signal out to cardiac_epochs
                outdata.data{carde_idx} = ensemble_add_data_struct_row(...
                  outdata.data{carde_idx},'subject_id',subid,'session',...
                  isess,'ensemble_id',sess.ensemble_id,'run',rnum,...
                  'trial',iep,'signal',{pksig});

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
              end
              
              % save final figure to file if not already saved
              if SAVE_EPOCHS && (m+1) ~= ep_per_fig
                print(pargs{:},figfname);
              end

              % save to outdata
              outdata.data{card_idx} = ensemble_add_data_struct_row(...
                  outdata.data{card_idx},'subject_id',lsubids,'session',...
                  lsessid,'ensemble_id',lesess,'run',lrun,'trial',...
                  [1:ns]','heart_rate',hr_bpm,'hr_variance',hr_var);
              
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
            if nrespd ~= ns
              warning(['data for response %s of different length (%d) '...
                  'than number of epochs (%d), these are not being '...
                  'added to the output\n'],rname,nrespd,ns);
              continue
            end

            % get indices for this channel and this response
            switch c
              case {'gsr','scr'}
                oidx = gsr_idx;
                ridx = gsr_cols.(rname);
              case {'cardiac','pulse'}
                oidx = card_idx;
                ridx = card_cols.(rname);
              otherwise
                warning('unknown channel type (%s)\n',c);
                continue
            end

            % save into outdata
            outdata.data{oidx}.data{ridx}(end-ns+1:end) = respdata;
            
          end % for iresp=1:length(resp_names
        end % for ic=1:nchan
      else
          
          warning(['no epoch start or end identified ... we could '...
              'treat the entire signal as one big epoch, but this '...
              'has not been coded yet']);
      end % if isstruct(epoch_start
    end % for irun=1:
  end % for isess = 1:

  if (GRAPH_EPOCHS && EPOCH_BY_RESP) || SAVE_DATA
    pfilt = struct();
    pfilt.include.all.subject_id = {subid};
    pfilt.include.all.path_type = {'physio_outdir'};
    lpathdata = ensemble_filter(pathdata,pfilt);
    if ~isempty(lpathdata.data{1})
      physdir = lpathdata.data{pcol.path}{1};
      check_dir(physdir);
    else
      warning(['could not find physio_outdir in pathdata, not graphing '...
          'epochs by response or saving data for subject %s'],subid);
      continue
    end

    % save data?
    if SAVE_DATA
      for ic=1:nchan
        c = channels{ic};
        % get channel idxs
        switch c
          case {'gsr','scr'}
            oidx = gsr_idx;
          case {'cardiac','pulse'}
            oidx = card_idx;
          otherwise
            warning('unknown channel type (%s)\n',c);
        end % switch c
        
        %%%% FIXME: temporary settings, remove before forgetting
        xsp.export = export;
        xsp.export.fname = fullfile(physdir,sprintf('%s_%s_export.txt',...
            subid,outdata.data{oidx}.type));
        xsp.sas = sas;
        xsp.sas.fname = fullfile(physdir,sprintf('%s_%s_export.sas',...
            subid,outdata.data{oidx}.type));
        xsp.sas.libname=sprintf('s%s_%s',subid,outdata.data{oidx}.type);
        ensemble_export_sastxt(outdata.data{oidx},xsp);
      end
    end % if SAVE_DATA

    if GRAPH_EPOCHS && EPOCH_BY_RESP
      figfname = fullfile(physdir,sprintf('%s_epochs_by_resp.ps',subid));
      % graph epochs by response type
      for iresp=1:length(resp_names)
        r = resp_names{iresp};
        for ic=1:nchan
          c = channels{ic};
          % get channel and response indexes
          switch c
            case {'gsr','scr'}
              oidx = gsr_idx;
              eidx = gsre_idx;
              ridx = gsr_cols.(rname);
              sidx = gsre_cols.signal;
            case {'cardiac','pulse'}
              oidx = card_idx;
              eidx = carde_idx;
              ridx = card_cols.(rname);
              sidx = carde_cols.signal;
            otherwise
              warning('unknown channel type (%s)\n',c);
              continue
          end % switch c

          % get mask matrix
          if iscell(outdata.data{oidx}.data{ridx}) && ...
                  all(cellfun(@isempty,outdata.data{oidx}.data{ridx}))
            warning('no responses for response_name %s\n',r);
            continue
          end

          [rm,ur] = make_mask_mtx(outdata.data{oidx}.data{ridx});
        
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
                sig=diff(pidxs);
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
