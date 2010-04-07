function ri = read_biopac(fname,params)

% Reads physiological data measured by the BioPac system, epochs data
% 
%   ri = read_biopac(fname,params)
% 
% The BioPac system has been used by the Janata Lab at the 3T
% scanner at the UC Davis Imaging Research Center, to collect physiological
% measures while participants were being scanned. This system, though,
% could conceivably be used to collect data outside of the scanner.
% Therefore, this script does not assume that physiological data it
% analyzes was collected in an NMR environment, and it does not assume that
% a TR is available for epoching purposes, but it does allow for a TR to be
% used for epoching, and it also allows for a more general presentation
% logfile-based epoching.
% 
% REQUIRES
%   - reading data -
%   fname - path to single file to be processed
%   params.channels(:) - struct array containing information about the
%     channels to be extracted and returned
%     .name - channel name, that will be converted to a fieldname in ri
%       i.e. if a channel name is "cardiac", then that data will be found
%       in ri.signal.cardiac in the output data
%     .comment_text - name given to the data stream in AcqKnowledge/BioPac,
%       as reflected in .hdr.per_chan_data.comment_text in the .acq file
%     NOTE: special channel name 'trigger' treated as scanner pulses,
%     assumed to be TTL pulses time-locked to volume acquisition onset
%   params.signal_bounds - 2-length vector containing start and end times
%       in seconds to initially bound the raw signal. this can be used in
%       the case that one physio file contains data from multiple runs.
%       Since data here are handled on a run-by-run basis, this can be used
%       to limit the initial data that is extracted to the samples
%       pertinent for the given run. a value of 0 in the first position
%       will be replaced with a 1, while a value of 0 in the second
%       position will be replaced with the final sample in the data set
%
%   - epoching -
%   params.bound_by_trigger - if this is true, and if there is a channel in
%     params.channels(:) with name 'trigger', and if this data is not empty
%     or a flat line, triggers will be used to define the epoch for this run
%     ... the epoch will be first trigger through last trigger+TR*srate ...
%     NOTE: if params.bound_by_trigger is set, and if a trigger channel
%     exists, params.link2pres will be ignored
%     NOTE: this requires params.TR ... if no params.TR is specified, it
%     will be estimated as median(diff(trigger_onsets_in_ms))/1000
%   params.baseline - the amount of time, in seconds, to allow at the
%     beginning of the signal for a baseline calculation.
%   params.TR - length (in seconds) between repetitions, for fMRI analyses.
%     see also: params.bound_by_trigger
%   params.link2pres - structure with settings to link up to a presentation 
%     response file ... MAKE SURE that you specify filt criteria that will 
%     return a set of items that then have a 1 to 1 correlation with event
%     onsets in a channel within the acq file. If the return from your
%     filtering criteria exactly matches the channel you've specified to 
%     link to, then all signals from the acq file will be aligned with the 
%     presentation data, and data falling before run_rel_time = 0 will be
%     thrown out
%     NOTE: if params.bound_by_trigger is set, and if a trigger channel
%     exists, link2pres will be ignored
%     NOTE: link2pres ASSUMES that the run begins AFTER physio data
%     collection begins. Even if presentation is started before physio data
%     collection, it will adjust all times in run_rel_time to what it
%     thinks is the beginning of the run, and all signals processed by
%     link2pres rely on this.
%     .presfname - path and filename of the presentation file to link to
%     .filt - filtering criteria applied to presfname data via
%       ensemble_filter ... you may want to use this struct to filter by
%       run, if you are calling read_biopac (or proc_physio_data) from
%       within a loop over runs, for instance
%     .channelname - ri.signal.(.channelname) will be compared to
%       filtered results from presfname. If a one-to-one match is found
%       between presfname filtered data and .channelname data, then epochs
%       will be defined by the beginning of the presentation run_rel_time,
%       and the last event returned by the presfname filtered data
%     .diff_tol - a diff will be done for all onsets given for presentation
%       data as well as link signal data, and then a diff will be done
%       between those diffs. If the diffs don't line up in length, or if
%       there is any diff between diffs greater than .diff_tol, link2pres
%       will no be applied
%     .bound_by_last_pulse - if set to 1, and if EVENT_TYPE={'Pulse'},
%       EVENT_CODE = {'255'} filtering criteria returns any rows, all data
%       after sample = ((lastPulseTime+TR*1000)/stime) will be dropped. If
%       no params.TR is specified, TR=median(diff(pulse_onsets_in_ms))/1000
%     .estimate_trigger - if a 'trigger' channel can not be found, but
%       presentation data is provided and contains events in the form of
%       EVENT_TYPE = {'Pulse'} && EVENT_CODE = 255, and if link2pres is
%       specified, after temporal linking is achieved, the pulse events
%       will be extracted and converted to triggers. Missing triggers will
%       be interpolated, either using the specified TR, or using a TR
%       estimated from the median pulse period.
% 
% OUTPUT
%   ri - 'runinfo' struct, in a form similar to that of read_keithley.m and
%     read_mate.m, containing fields for data streams, events, and onsets,
%     also containing a .meta field, with .meta.srate (sampling rate),
%     .meta.stime (sampling time, or # ms per sample), and 
%     .meta.baseline_samps (# of baseline samples). This differs from that
%     of read_keithley and read_mate in that the actual signal fields are
%     found in ri.signal, rather than within ri
%       NOTE: ri.trigger values are sample indices, NOT times in ms or sec
% 
% FB 2008.10.12
% FB 2009.06.24 - added baseline calculation

ri = struct('signal',struct(),'meta',struct());
global r

r = init_results_struct;
r.type = 'physio_biopac';  % Identify the type of this reporting instance
r.report_on_fly = 1;

if ~isempty(fname) && exist(fname,'file')
  % read data
  pdata = load_acq(fname);
  chans = {pdata.hdr.per_chan_data(:).comment_text};
else
  % data can't be found, skip this file
  msg = sprintf('couldn''t find biopac file %s, SKIPPING\n',fname);
  r = update_report(r,msg);
  return
end

% get sampling rate
stime = pdata.hdr.graph.sample_time; % msec per sample
srate = 1000/stime; % 1000 msec / stime = sampling rate
ri.meta.srate = srate;
ri.meta.stime = stime;

% get TR
if isfield(params,'TR')
  TR = params.TR;
else
  TR = 0;
  warning('no TR found in the params struct\n');
end

% signal bounds?
if isfield(params,'signal_bounds') && length(params.signal_bounds) == 2
  % convert from params.signal_bounds in secs to sig_start/end in samples
  sig_start = params.signal_bounds(1)*srate;
  sig_end = params.signal_bounds(2)*srate;
  if ~isnumeric(sig_start) || sig_start < 1
    sig_start = 1;
  end
  if ~isnumeric(sig_end) || sig_end < 1 || sig_end > size(pdata.data,1)
    sig_end = size(pdata.data,1);
  end
else
  sig_start = 1;
  sig_end = size(pdata.data,1);
end

% find channels, get data
ctxts = {pdata.hdr.per_chan_data.comment_text};
if isfield(params,'channels')
  nc = length(params.channels);
  for ic=1:nc
    ctxt = params.channels(ic).comment_text;
    name = params.channels(ic).name;
    cidx = strmatch(ctxt,ctxts);
    if ~isempty(cidx)
      ri.signal.(name) = pdata.data(sig_start:sig_end,cidx);
    end
  end
end

% find trigger onsets, if they exist
if isfield(ri.signal,'trigger')
  % onsets and offsets are in samples, NOT time in msec or sec
  [ri.trigger.onsets,ri.trigger.offsets] = ...
      find_thresh_cross(ri.signal.trigger,params);
end

% extract run-length epochs
% epoch_start and epoch_end are expressed in samples, NOT in time
if isfield(ri,'trigger') && isfield(ri.trigger,'onsets') ...
      && ~isempty(ri.trigger.onsets) && isfield(params,'bound_by_trigger') ...
      && params.bound_by_trigger
  % check for TR
  if ~exist('TR','var') || ~TR
    warning('estimating TR from trigger onsets');
    TR = median(diff(ri.trigger.onsets*stime))/1000;
    warning('TR estimated to be %ds',TR);
  end

  % bound by triggers
  epoch_start = ri.trigger.onsets(1);
  epoch_end = ri.trigger.onsets(end)+TR*srate; % add one TR at the end

elseif isfield(params,'link2pres') ...
        && isfield(params.link2pres,'channelname')...
        && ischar(params.link2pres.channelname) ...
        && isfield(ri.signal,params.link2pres.channelname)
    
  % extract signal
  linksig = ri.signal.(params.link2pres.channelname);
  [lsons,lsoffs] = find_thresh_cross(linksig,params);
  lsdiff_ms = diff(lsons)*stime;
  nlsdiff = length(lsdiff_ms);
  
  % link to presentation file, bound by presentation data
  pfname = params.link2pres.presfname;
  presdata = load(pfname);
  pcol = set_var_col_const(presdata.vars);
  
  if isfield(params.link2pres,'filt')
    pfilt = params.link2pres.filt;
    linkdata = ensemble_filter(presdata,pfilt);
  else
    linkdata = presdata;
  end
  pdiff_ms = diff(linkdata.data{pcol.RUN_REL_TIME});
  npdiff = length(pdiff_ms);

  if isfield(params.link2pres,'diff_tol')
    diff_tol = params.link2pres.diff_tol;
  else
    diff_tol = 100; % 100ms tolerance
  end
  nlsd = length(lsdiff_ms);
  npd  = length(pdiff_ms);

  % the first two conditionals within this if/else train are a hack to
  % account for an additional leading or trailing signal in either linksig
  % or presdata, that doesn't appear in the other ... essentially, if
  % linksig has items 1 2 3 and pdata has 2 3 or 1 2 (i.e. pdata is missing
  % a leading or trailing element, but otherwise is intact), or vice versa
  % (linksig is missing an element), we can still link the two measurements
  % together
  
%   FIXME: right now, d3 = d1>d2, and the sum of # of indices where one is
%   greater than the other is used to figure out where 'less difference'
%   exists, but this is not really a metric of 'less difference', rather
%   'indices with less difference' ... might another formulation work
%   better? how about d4 = d1-d2, then compare abs(sum(d4(d3))) to
%   abs(sum(d4(~de))) ... ? if problems arise, I will try this strategy
  if (nlsd - npd) == 1
    % one more event in the link signal than in the presentation file
    
    % assume the first pres signal is missing
    d1 = abs(diff([lsdiff_ms'; [lsdiff_ms(1); pdiff_ms]']));
    
    % assume the last pres signal is missing
    d2 = abs(diff([lsdiff_ms'; [pdiff_ms; lsdiff_ms(end)]']));
    
    % greater diff for first signal compared to last signal assumption?
    d3 = d1 > d2;
    
    if (sum(d3) < nlsd/2)
      % smaller difference between link sig and pres times when assuming
      % the first pres time is missing
      % in other words, linksig(2) = prestime(1)
      if ~any(d1 > diff_tol)
        % no differences greater than diff_tol, safe to assume the first
        % pres time is missing
        pbegin = linkdata.data{pcol.RUN_REL_TIME}(1); % 1st pres time        
        % check to see that the first epoch time is not before the
        % beginning of the run as presentation has calculated it
        if lsons(2)*stime - lsons(1)*stime > pbegin
          msg = sprintf(['possible error calculating the beginning of '...
              'the epoch. There are more link signal events %d than '...
              'presentation linking events %d, and it appears that the '...
              'first presentation event is missing, but if this is the '...
              'case, the first link signal event is before the start '...
              'of the run, when comparing event times to presentation '...
              'run-relative times'],nlsd,npd);
          r = update_report(r,msg);
          ri = struct();
          return
        else
          % start of the epoch is the onset of the second link signal
          % event, minus the time of the first (but really, the second)
          % presentation linking event, converted to link signal samples
          epoch_start = round((lsons(2)*stime - pbegin)/stime);
        end
      else
        % differences greater than diff_tol, likely a middle time that is
        % missing
        msg = sprintf('missing presentation link, can not compute');
        r = update_report(r,msg);
        ri = struct();
        return
      end
    elseif (sum(d3) > nlsd/2)
      % smaller difference between link sig and pres times when assuming
      % the last pres time is missing
      % in other words, linksig(1) = prestime(1)
      if ~any(d2 > diff_tol)
        % no differences greater than diff_tol, safe to assume the last
        % pres time is missing
        pbegin = linkdata.data{pcol.RUN_REL_TIME}(1); % 1st pres time
        % start of the epoch is the onset of the first link signal event,
        % minus the time of the first presentation linking event, converted
        % to link signal samples
        epoch_start = round((lsons(1)*stime - pbegin)/stime);
      else
        % differences greater than diff_tol, likely a middle time that is
        % missing
        msg = sprintf('missing presentation link, can not compute');
        r = update_report(r,msg);
        ri = struct();
        return
      end
    else
      msg = sprintf('missing presentation link, can not compute');
      r = update_report(r,msg);
      ri = struct();
      return
    end % if (sum(d3) < nlsd/2
  elseif (npd - nlsd) == 1
    % one more event in the presentation file than in the link signal
    
    % assume the first link signal is missing
    d1 = abs(diff([pdiff_ms'; [pdiff_ms(1); lsdiff_ms]']));
    
    % assume the last link signal is missing
    d2 = abs(diff([pdiff_ms'; [lsdiff_ms; pdiff_ms(end)]']));
    
    % greater diff for first signal compared to last signal assumption?
    d3 = d1 > d2;
    
    if (sum(d3) < npd/2)
      % smaller difference between link sig and pres times when assuming
      % the first link signal is missing
      % in other words, prestime(2) == linksig(1)
      if ~any(d1 > diff_tol)
        % no differences greater than diff_tol, safe to assume the first
        % link signal is missing
        pbegin = linkdata.data{pcol.RUN_REL_TIME}(2); % 2nd pres time
        % check to see that the beginning of physio data collection is
        % before the beginning of the run as presentation has logged it
        if lsons(1)*stime < pbegin
          msg = sprintf(['possible error calculating the beginning of '...
              'the epoch. There are more presentation linking events %d'...
              ' than link signal events %d, and it appears that the '...
              'first link signal event is missing, but if this is the '...
              'case, the first presentation linkin event is before the '...
              'start of physio data collection, when comparing physio '...
              'data to presentation run-relative times'],npd,nlsd);
          r = update_report(r,msg);
          ri = struct();
          return
        else
          % start of the epoch is the onset of the first link signal
          % event, minus the time of the second presentation linking event,
          % converted to link signal samples
          epoch_start = round((lsons(1)*stime - pbegin)/stime);
        end
      else
        % differences greater than diff_tol, likely a middle time that is
        % missing
        msg = sprintf('missing link signal, can not compute');
        r = update_report(r,msg);
        ri = struct();
        return
      end
    elseif (sum(d3) > npd/2)
      % smaller difference between link sig and pres times when assuming
      % the last link signal is missing
      % in other words, linksig(1) = prestime(1)
      if ~any(d2 > diff_tol)
        % no differences greater than diff_tol, safe to assume the last
        % link signal is missing
        pbegin = linkdata.data{pcol.RUN_REL_TIME}(1); % 1st pres time
        % start of the epoch is the onset of the first link signal event,
        % minus the time of the first presentation linking event, converted
        % to link signal samples
        epoch_start = round((lsons(1)*stime - pbegin)/stime);
      else
        % differences greater than diff_tol, likely a middle time that is
        % missing
        msg = sprintf('missing link signal, can not compute');
        r = update_report(r,msg);
        ri = struct();
        return
      end
    else
      msg = sprintf('missing link signal, can not compute');
      r = update_report(r,msg);
      ri = struct();
      return
    end
  elseif npd ~= nlsd
    msg = sprintf(['link signal size (%d) is not equal in length to '...
        'presentation data size (%d), therefore link2pres skipped\n'],...
        nlsdiff,npdiff);
    r = update_report(r,msg);
    ri = struct();
    return
  elseif any(diff([lsdiff_ms'; pdiff_ms']) > diff_tol)
    msg = sprintf(['more than %d ms diff btwn link signal onset diffs '...
        'and presentation data onset diffs, therefore link2pres skipped\n'],...
        diff_tol);
    r = update_report(r,msg);
    ri = struct();
    return
  else
    pbegin = linkdata.data{pcol.RUN_REL_TIME}(1);
    epoch_start = round((lsons(1)*stime - pbegin)/stime);
  end

  epoch_end = length(linksig);
  
  if isfield(params.link2pres,'bound_by_last_pulse') || ...
          isfield(params.link2pres,'estimate_trigger')
    lpfilt.include.all.EVENT_TYPE={'Pulse'};
    lpfilt.include.all.EVENT_CODE={'255'};
    lpfilt.include.all.RUN=params.run;
    lpdata = ensemble_filter(presdata,lpfilt);
    rrtimes = lpdata.data{pcol.RUN_REL_TIME};
    
    if ~isempty(lpdata.data{1})
      if ~exist('TR','var') || ~TR
        warning('estimating TR from pulse data');
        TR = median(diff(rrtimes))/1000;
        warning('TR estimated to be %ds',TR);
      end
      
      if isfield(params.link2pres,'estimate_trigger')
        cpp = check_pulse_periods(rrtimes);
        pptl = find(cpp.data.pulse_period_too_long_mask);
        for ipp=1:length(pptl)
          intpp = rrtimes(pptl(ipp)):(TR*1000):rrtimes(pptl(ipp)+1);
          rrtimes = union(rrtimes,intpp);
        end
        samples = rrtimes/stime;
        ri.trigger.onsets = samples;
        ri.trigger.offsets = samples+1;
      end
      
      if isfield(params.link2pres,'bound_by_last_pulse')
        % the end of the epoch = run_rel_time/stime + epoch_start
        pend = ceil((lpdata.data{pcol.RUN_REL_TIME}(end)+TR*1000)/stime);
        epoch_end = pend + epoch_start;
      end
    else
      warning('no pulses found, no trigger estimation or bounding');
    end
  end

else
  warning('could not generate run-length epochs');
  return
end

% add baseline period?
if isfield(params,'baseline')
    blstart = epoch_start - params.baseline*srate;
    if blstart < 1, blstart = 1; end
    ri.meta.baseline_samps = epoch_start - blstart;
    epoch_start = blstart;
else
    ri.meta.baseline_samps = 0;
end

% apply bounds
flds = fieldnames(ri.signal);
nflds = length(flds);
for ifld = 1:nflds
    fld = flds{ifld};
    nsamp = length(ri.signal.(fld));
    if epoch_start > nsamp
      continue
    end
    if epoch_end > nsamp
      epoch_end = nsamp;
    end
    ri.signal.(fld) = ri.signal.(fld)(epoch_start:epoch_end);
end

% adjust triggers?
if isfield(ri,'trigger') && isfield(ri.trigger,'onsets') ...
      && ~isempty(ri.trigger.onsets) && isfield(params,'bound_by_trigger') ...
      && params.bound_by_trigger
  ri.trigger.onsets = ri.trigger.onsets + ri.meta.baseline_samps;
  ri.trigger.offsets = ri.trigger.offsets + ri.meta.baseline_samps;
end
