function outdata = ensemble_fmri_physio_filter(indata,defs)

% filters pre-processed physiological data, stores .mat file on disk, 1 per run
% 
% REQUIRES
%   defs.physio_filter.target_channels
%   defs.physio_filter.SIG_CHECK
%   defs.physio_filter.butter.order
%   defs.physio_filter.butter.cutoff
%   defs.figs.print_args
%   defs.fmri
% 
% 2008/10/06 FB - started coding

outdata = ensemble_init_data_struct();

global r

r = init_results_struct;

r.type = 'fmri_physio';  % Identify the type of this reporting instance
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
        paths = indata{idata};
        pscol = set_var_col_const(paths.vars);
      case {'physio_paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
      case {'sinfo'}
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
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
  if exist('paths','var') && ~isempty(paths.data{1})
    if length(nsub_proc) == 1
      pfilt = struct();
      pfilt.include.all.subject_id = proc_subs;
      lpathdata = ensemble_filter(paths,pfilt);
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

outdata.vars = [outdata.vars 'physio_paths'];
ppaths_idx = length(outdata.vars);
outdata.data{ppaths_idx} = ensemble_init_data_struct();
outdata.data{ppaths_idx}.type='physio_paths';
outdata.data{ppaths_idx}.vars = {'subject_id','session',...
    'ensemble_id','run','path'};
outdata.data{ppaths_idx}.data{1} = {};
outdata.data{ppaths_idx}.data{2} = [];
outdata.data{ppaths_idx}.data{3} = [];
outdata.data{ppaths_idx}.data{4} = [];
outdata.data{ppaths_idx}.data{5} = {};
ppcol = set_var_col_const(outdata.data{ppaths_idx}.vars);

target_channels = defs.physio_filter.target_channels;
try SIG_CHECK = defs.physio_filter.SIG_CHECK; catch SIG_CHECK = 0; end
try b_order = defs.physio_filter.butter.order; catch b_order = 0; end
try b_order = defs.physio_filter.butter.order; catch b_order = 12; end
try b_cutoff = defs.physio_filter.butter.cutoff; catch b_cutoff = 0.5; end

if ~isfield(defs,'figs') || ~isfield(defs.figs,'printargs') || ...
        ~iscell(defs.figs.printargs) || isempty(defs.figs.printargs)
  defs.figs.printargs = {'-dpsc','-append'};
end

if b_order
  [fb,fa] = butter(b_order,b_cutoff);
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
    
    session_stub = sess.id;
    r = update_report(r,sprintf('\t\t\t%s\n', session_stub));
    
    % init session vars
    physio = sinfo(isub).sessinfo(isess).physio;
    pparams = init_physmon_params(physio.source);
    pparams.channels = physio.channels;
    if isfield(physio,'link2pres')
      pparams.link2pres = physio.link2pres;
    end
    
    % Obtain information about the scanning protocol that was used in this
    % session
    protocol_idx = strmatch(sess.protocol_id, {defs.fmri.protocol.id}, ...
	'exact');
    
    pparams.TR = defs.fmri.protocol(protocol_idx).epi.tr;
    pparams.nslice_per_vol = defs.fmri.protocol(protocol_idx).epi.nslices;
    
    % get physio_indir/outdir, behav_outdir
    sfilt = struct();
    sfilt.include.all.subject_id = {subid};
    sfilt.include.all.session = isess;
    spaths = ensemble_filter(pathdata,sfilt);
    
    % Determine how many runs we're dealing with
    runs = sess.use_epi_runs;
    nruns = length(runs);  
    
    %
    % START OF RUN LOOP
    %
    for irun = 1:nruns

      rnum = runs(irun);

      msg = sprintf('\t\tPROCESSING RUN %d (%d/%d)\n',rnum,irun,nruns);
      r = update_report(r,msg);
        
      rfilt = sfilt;
      rfilt.include.all.run = rnum;
      rpaths = ensemble_filter(spaths,rfilt);
      
      % get physio data
      physiofname = rpaths.data{pcol.path}{1};
      if isempty(physiofname)
        msg = sprintf('no physiofname found, skipping run %d (%d/%d)\n',...
            rnum,irun,nruns);
        r = update_report(r,msg);
        continue
      end
      
      [fp,fn,fx] = fileparts(physiofname);
      EEG = pop_loadset('filename',sprintf('%s%s',fn,fx),'filepath',fp);
      tgt_chan_idxs = find(ismember({EEG.chanlocs(:).labels},target_channels));
      nchans = length(tgt_chan_idxs);

      EEGraw = EEG; % save a copy of the original dataset
      
      for ic=1:nchans
        tci = tgt_chan_idxs(ic);
        cname = EEG.chanlocs(tci).labels;
        fprintf(1,'chan %d (%d), ',ic,tci);

        % get data for this channel
        Y = EEG.data(tci,:);

        if ~isempty(strmatch(cname,{'scr','respir'}))
        
          %
          % remove any scanner-related artifact in scr and respir channels
          % by calculating mean wave-form for each TR, and regressing the
          % channel signal on a design matrix based on this mean waveform
          %

          % create TR event channel
          fprintf(1,'creating TR event channel\n');
          tr_idxs = 1:pparams.TR*EEG.srate:EEG.pnts; % in samples
          ntr = length(tr_idxs);
          tr_events = [ones(ntr,1) (tr_idxs/EEG.srate)']; % in seconds
          EEGe = pop_importevent(EEG,'append','no','event',tr_events,...
              'fields',{'type','latency'},'timeunit',1);

          % extract TR epochs
          fprintf(1,'extracting epochs\n');
          EEGe = pop_epoch(EEGe,{},[0 pparams.TR],'epochinfo','yes');
          E = EEGe.data(tci,:);

          % calculate mean TR waveform for each channel
          fprintf(1,'calculating mean TR waveform for %s channel, ',cname);
          W = [];
          for it=1:EEGe.trials
            W = [W; E(:,it)];
          end
          mW = mean(W); % mean waveform across TRs

          % create regression design matrix
          dM = zeros(EEG.pnts,EEGe.trials);
          fprintf(1,'dM, ');
          for it=1:EEGe.trials
            start = (EEGe.pnts*(it-1)+1);
            stop  = EEGe.pnts*it;
            dM(start:stop,it) = mW;
          end
          
          % are there any samples not accounted for at the end of dM?
          if stop < size(dM,1)
            if (size(dM,1) - (stop - 1)) > length(mE)
              % left-over size is greater than the length of a TR
              error('there are less trials than there should be\n');
            else
              % fill out extra time at the end of the design matrix
              start = stop + 1;
              stop = size(dM,1);
              dM(start:stop,end+1) = repmat(mE(1:stop-start),[it 1]);
            end
          end
          dM(:,end+1) = 1; % must place a constant as the last regressor

          % filter using filtfilt?
          if b_order
            fprintf(1,'filtering, ');
            Y = filtfilt(fb,fa,Y);
          end
          
          fprintf(1,'regressing, ');
          [b,bint,r,rint,stats] = regress(EEGraw.data,dM(:,:));
          Y = Y - r';
          
          if SIG_CHECK
            logfname = fullfile(fp,sprintf('%s_sigcheck.ps',fn));
            figure();

            subplot(3,1,1);
            plot(EEG.data(tci,:));
            title(sprintf('raw %s',cname));

            subplot(3,1,2);
            plot(r);
            title(sprintf('residual %s',cname));

            subplot(3,1,3);
            plot(Y);
            title(sprintf('processed (raw - residual) %s',cname));

            print(logfname,defs.figs.printargs{:});
          end
        else
          % filter using filtfilt?
          if b_order
            fprintf(1,'filtering, ');
            Y = filtfilt(fb,fa,Y);
          end
          
          if SIG_CHECK
            logfname = fullfile(fp,sprintf('%s_sigcheck.ps',fn));
            figure();

            subplot(2,1,1);
            plot(EEG.data(tci,:));
            title(sprintf('raw %s',cname));

            subplot(2,1,2);
            plot(Y);
            title(sprintf('processed %s',cname));
            print(logfname,defs.figs.printargs{:});
          end                    
        end

        EEGraw.data(tci,:) = Y;
      end % for ic=1:nchans

      set_fname = sprintf('%s_filtered%s',fn,fx);
      pop_saveset(EEGraw,'filename',set_fname,'filepath',fp);
        
      set_fpn = fullfile(fp,set_fname);

      outdata.data{ppaths_idx}.data{ppcol.subject_id} = ...
          [outdata.data{ppaths_idx}.data{ppcol.subject_id}; subid];
      outdata.data{ppaths_idx}.data{ppcol.session} = ...
          [outdata.data{ppaths_idx}.data{ppcol.session}; isess];
      outdata.data{ppaths_idx}.data{ppcol.ensemble_id} = ...
          [outdata.data{ppaths_idx}.data{ppcol.ensemble_id}; sess.ensemble_id];
      outdata.data{ppaths_idx}.data{ppcol.run} = ...
          [outdata.data{ppaths_idx}.data{ppcol.run}; irun];
      outdata.data{ppaths_idx}.data{ppcol.path} = ...
          [outdata.data{ppaths_idx}.data{ppcol.path}; set_fpn];
    end % for irun
    
  end % for isess
end % for isub=