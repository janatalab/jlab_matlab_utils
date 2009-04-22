function outdata = ensemble_fmri_physio_summ(indata,defs)

% calculates summary statistics for physio data collected during fMRI
% 
% FB 2009.04.21

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'fmri_physio_summ';  % Identify the type of this reporting instance
r.report_on_fly = 1;

%%%% FIXME: allow for a param?
% load colormap
load new_seismic;
colormap(new_seismic);

% init EEGLAB
eeglab('initpaths');

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

outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

if isfield(defs,'physio_summ')
  ps = defs.physio_summ;
else
  error('please provide params through defs.physio_summ\n');
end

try channels = ps.channels; catch error('must provide channel names\n'); end
try int_begin = ps.spc_integral.begin; catch int_begin = 2000; end
try int_end   = ps.spc_integral.end; catch int_end = 6000; end
try pk.thresh = ps.peaks.thresh; catch pk.thresh = 0; end
try pk.calcfrom = ps.peaks.calcFrom; catch pk.calcfrom = 'left'; end
try pkh.thresh = ps.peak_heights.thresh; catch pkh.thresh = 0.75; end
try proxlim = ps.proxlim; catch proxlim = 0.6; end
try GRAPH_EPOCHS = ps.GRAPH_EPOCHS; catch GRAPH_EPOCHS = 0; end
try SC_CONTOUR = ps.SC_CONTOUR; catch SC_CONTOUR = 0; end
try export = ps.export; catch export.write2file = 0; end
try sas = ps.sas; catch sas.libsave = 0; end
try epoch_start = ps.epoch_start; catch epoch_start =0; end
try epoch_end = ps.epoch_end; catch epoch_end = 0; end
try responses = ps.responses; catch responses = 0; end

if isstruct(responses)
  resp_names = {responses.name};
else
  resp_names = {};
end

supported = {'gsr','cardiac'};
nchan = length(channels);
nmatch = 0;

for ic=nchan:1
  if ~isempty(strmatch(channels{ic},supported))
    match = match + 1;
  end
end

if ~nmatch, error('no supported channels provided'), end

eeglab('initpaths');

if ~isempty(strmatch('gsr',channels)) || ~isempty(strmatch('scr',channels))
    % initialize a gsr summary data output structure
    outdata.vars = [outdata.vars 'gsr_summ'];
    gsr_idx = length(outdata.vars);
    outdata.data{gsr_idx} = ensemble_init_data_struct();
    outdata.data{gsr_idx}.type = 'gsr_summ';
    outdata.data{gsr_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','integral','mean_activation',...
        'ttl_integral','ttl_mean_activation','npeaks','first_peak_time',...
        'first_peak_height',resp_names{:}};
    gsr_cols = set_var_col_const(outdata.data{gsr_cols}.vars);
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
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_time} = [];
    outdata.data{gsr_idx}.data{gsr_cols.first_peak_height} = [];
    for in=1:length(resp_names)
      outdata.data{gsr_idx}.data{gsr_cols.(resp_names{i})} = [];
    end
end

if ~isempty(strmatch('cardiac',channels))
    % initialize a cardiac summary data output structure
    outdata.vars = [outdata.vars 'cardiac_summ'];
    card_idx = length(outdata.vars);
    outdata.data{card_idx} = ensemble_init_data_struct();
    outdata.data{card_idx}.type = 'cardiac_summ';
    outdata.data{card_idx}.vars = {'subject_id','session','ensemble_id',...
        'run','trial','npeaks','heart_rate','hr_variance',resp_names{:}};
    card_cols = set_var_col_const(outdata.data{card_cols}.vars);
    outdata.data{card_idx}.data{card_cols.subject_id} = {};
    outdata.data{card_idx}.data{card_cols.session} = [];
    outdata.data{card_idx}.data{card_cols.ensemble_id} = [];
    outdata.data{card_idx}.data{card_cols.run} = [];
    outdata.data{card_idx}.data{card_cols.trial} = [];
    outdata.data{card_idx}.data{card_cols.npeaks} = [];
    outdata.data{card_idx}.data{card_cols.heart_rate} = [];
    outdata.data{card_idx}.data{card_cols.hr_variance} = [];
    for in=1:length(resp_names)
      outdata.data{card_idx}.data{card_cols.(resp_names{i})} = [];
    end
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
            ~exist(pdata.data{prescol.path}{1})
      error('can not find presentation data for subject %s, session %d',...
          subid,isess)l
    else
      presdata = load(pdata.data{prescol.path}{1});
      presdcol = set_var_col_const(presdata.vars);
    end
    
    nruns = length(sess.use_epi_runs);
    for irun=1:nruns
      
      msg = sprintf('\t\tPROCESSING RUN (%d/%d)\n', irun,nruns);
      r = update_report(r,msg);
      
      pfilt.include.all.subject_id = {subid};
      pfilt.include.all.session = isess;
      pfilt.include.all.run = sess.use_epi_runs(irun);
      
      lphys = ensemble_filter(physio_paths,pfilt);
      if isempty(lphys.data{physcol.path}) || ...
              ~exist(lphys.data{physcol.path}{1},'file')
        warning('no physio file found for %s, sess %d, run %d\n',...
            subid,isess,sess.use_epi_runs(irun));
        continue
      end
      
      [fp,fn,fx] = fileparts(lphys.data{physcol.path}{1});
      EEG = pop_loadset('filename',sprintf('%s%s',fn,fx),'filepath',fp);
      
      % epoch the EEG data, set -200ms or -250ms start
      % use epoch_start and epoch_end, get from presdata, line up starts
      % and ends
      
      % extract each response type, line up with starts/ends
      
      % iterate over channels
      % calculate summary stats for each channel
      % save to outdata
      
      %%%% for heart rate
      pidxs = find_peaks(EEG.data(ichan,:,itrial),pk); %% or is it pkh?
      diffs = pidxs(2:end) - pidxs(1:end-1);
      diffs(diffs < 4*std(diffs)) = []; % remove incredibly short diffs
      heart_rate_Hz = 1/(mean(diffs)/EEG.srate); % heart rate in Hz
      heart_rate_bpm = 1/(mean(diffs)/EEG.srate) * 60; % heart rate, bpm
      hr_variance = var(diffs);
      npeaks = length(pidxs);
    end % for irun=1:
  end % for isess = 1:
end % for isub=1:


