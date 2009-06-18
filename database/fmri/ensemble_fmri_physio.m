function outdata = ensemble_fmri_physio(indata,defs)

% processes physiological data, stores .mat file on disk, 1 per run
% 
% If you are linking physio data to presentation file data, make sure to
% process the presentation file data first, so that it is available to this
% script.
% 
% This script does not preprocess or normalize, but it does cut up and
% limit per run.
% 
% NOTE: this script also exports data to EEGLAB .set format
% NOTE: this script assumes a slightly different format for the 'ri' physio
% data format that has been used in the past, in read_keithley and
% read_mate ... whereas in those scripts, channels were saved (such as
% 'cardiac' and 'respir') as fields of the main 'ri' struct, but in the new
% 'modified' format, the data fields are saved in ri.signal
% (ri.signal.respir, ri.signal.cardiac), so that meta data and header
% information can be saved as fields of 'ri'
% 
% REQUIRES
%   defs.sinfo(:).sessinfo(:).physio
%   defs.sinfo(:).sessinfo(:).physio.source
%   defs.sinfo(:).sessinfo(:).physio.resp_titles
%   defs.physio.SIG_CHECK
%   defs.physio.LINK2FMRI
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

outdata.vars = [outdata.vars 'physio_paths'];
ppaths_idx = length(outdata.vars);
outdata.data{ppaths_idx} = ensemble_init_data_struct();
outdata.data{ppaths_idx}.type='physio_paths';
outdata.data{ppaths_idx}.vars = {'subject_id','session',...
    'ensemble_id','run','path','preprocessed'};
outdata.data{ppaths_idx}.data{1} = {};
outdata.data{ppaths_idx}.data{2} = [];
outdata.data{ppaths_idx}.data{3} = [];
outdata.data{ppaths_idx}.data{4} = [];
outdata.data{ppaths_idx}.data{5} = {};
outdata.data{ppaths_idx}.data{6} = [];
ppcol = set_var_col_const(outdata.data{ppaths_idx}.vars);

ecdparams.outDataName = 'physio_data';

try SIG_CHECK = defs.SIG_CHECK; catch SIG_CHECK = 1; end
try LINK2FMRI = defs.LINK2FMRI; catch LINK2FMRI = 0; end

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
    
    if LINK2FMRI
      % Obtain information about the scanning protocol that was used in this
      % session
      protocol_idx = strmatch(sess.protocol_id, {defs.fmri.protocol.id}, ...
          'exact');
      pparams.TR = defs.fmri.protocol(protocol_idx).epi.tr;
      pparams.nslice_per_vol = defs.fmri.protocol(protocol_idx).epi.nslices;
    end
    
    
    % get physio_indir/outdir, behav_outdir
    sfilt = struct();
    sfilt.include.all.subject_id = {subid};
    sfilt.include.all.session = isess;
    spaths = ensemble_filter(pathdata,sfilt);
    
    indx = strmatch('physio_indir',spaths.data{pcol.path_type});
    physio_indir = spaths.data{pcol.path}{indx};

    outdx = strmatch('physio_outdir',spaths.data{pcol.path_type});
    physio_outdir = spaths.data{pcol.path}{outdx};
    
    if isfield(physio,'link2pres')
      behdx = strmatch('behav_outdir',spaths.data{pcol.path_type});
      behav_outdir = spaths.data{pcol.path}{behdx};
    end
    
    % Determine how many runs we're dealing with
    if LINK2FMRI
      runs = sess.use_epi_runs;
    else
      if isfield(sess,'use_physio_runs')
        runs = sess.use_physio_runs;
      else
        runs = 1:size(physio.datafiles,1);
      end
    end
    nruns = length(runs);  
    
    %
    % START OF RUN LOOP
    %
    for irun = 1:nruns

      msg = sprintf('\t\tPROCESSING RUN (%d/%d)\n', irun,nruns);
      r = update_report(r,msg);
        
      % get physio data
      physiofname = physio.datafiles{runs(irun),1};
      if isempty(physiofname)
        msg = sprintf('no physiofname found, skipping run %d\n',irun);
        r = update_report(r,msg);
        continue
      end
      physiopath = fullfile(physio_indir,physiofname);
      if ~exist(physiopath,'file')
        msg = sprintf('physiopath doesn''t exist, skipping run %d\n',irun);
        r = update_report(r,msg);
        continue
      end
      lparams = pparams;
      rparams = physio.datafiles{runs(irun),2};
      if isfield(rparams,'link2pres')
        lparams.run=irun;
        lparams.link2pres = rparams.link2pres;
        lparams.link2pres.presfname = fullfile(behav_outdir,...
            sprintf('%s_sess%d_present.mat',subid,isess));
        lparams.link2pres.filt.include.any.RUN = irun;
      end
      if isfield(rparams,'signal_bounds')
        lparams.signal_bounds = rparams.signal_bounds;
      end
      
      % read-in physio data
      ri = proc_physio_data(physiopath,lparams);

      if isempty(fields(ri))
        continue
      end
      
      % import into EEGLAB
      snames = fields(ri.signal);
      nchans = length(snames);
      npts   = length(ri.signal.(snames{1}));
      signal = zeros(nchans,npts);
      chanlocs = struct();
      for ic=1:nchans
        signal(ic,:) = ri.signal.(snames{ic});
        chanlocs(ic).labels = snames{ic};
      end
      
      EEG = pop_importdata('dataformat','array','data',signal,'srate',...
          ri.meta.srate,'subject',subid,'session',sess.ensemble_id,...
          'chanlocs',chanlocs,'nbchan',nchans,'pnts',npts);
      
      set_fname = sprintf('%s_sess%d_run%d_physio.set',subid,isess,irun);
      pop_saveset(EEG,'filename',set_fname,'filepath',physio_outdir);
      set_fpn = fullfile(physio_outdir,set_fname);

      outdata.data{ppaths_idx} = ensemble_add_data_struct_row(...
          outdata.data{ppaths_idx},'subject_id',subid,'session',isess,...
          'ensemble_id',sess.ensemble_id,'run',irun,'path',set_fpn,...
          'preprocessed',0);
      
      if SIG_CHECK
        sigs = {EEG.chanlocs(:).labels};
        nsig = length(sigs);
        figure();
        for i=1:nsig
          sig = sigs{i};
          subplot(nsig,1,i);
          plot(EEG.data(i,:));
          ylabel(sig);
        end % for i=1:nsig
        if isfield(defs,'figs') && isfield(defs.figs,'write2file') && ...
                defs.figs.write2file
            figfname = fullfile(physio_outdir,sprintf(...
                '%s_sess%d_run%d_physioSigCheck.ps',subid,isess,irun));
            print(figfname,defs.figs.printargs{:});
        end
      end % if SIG_CHECK
      
      %%%% Filter this run data?
      if ~isempty(strmatch('filter',fieldnames(lparams.channels)))
        EEGf = EEG;
        for ic=1:length(lparams.channels)
          if isfield(lparams.channels(ic),'filter')
            % get EEG channel index for this channel
            cidx = strmatch(lparams.channels(ic).name,...
                {EEGf.chanlocs(:).labels});
            if length(cidx) > 1
              cidx = cidx(1);
            elseif isempty(cidx)
              warning(['no match between lparams.channel %d and EEG ',...
                  'struct!? skipping\n'],ic);
              continue
            end % if length(cidx
          
            % filter data?
            lf = lparams.channels(ic).filter;
            if ~isstruct(lf)
              warning(['filter struct for channel %d is not a struct, '...
                  'skipping\n'],ic);
            else
              for icf=1:length(lf)
                if ~isfield(lf(icf),'fun')
                  warning(['filter struct %d for channel %d does not '...
                      'contain a function handle (.fun field), '...
                      'skipping\n'],icf,ic);
                else
                  lfh = parse_fh(lf(icf).fun);
                  fname = func2str(lfh);

                  cparams = lparams;
                  cparams.(fname) = lf(icf);
                  EEGf = lfh(EEGf,cparams);
                end % if ~isfield(lf(icf
              end % for icf=1:
            end % if ~isstruct(lf
          end % if isfield(lparams.channels(ic
        end % for ic=1:length(lparams.ch

        set_fname = sprintf('%s_sess%d_run%d_physio_filtered.set',subid,isess,irun);
        pop_saveset(EEGf,'filename',set_fname,'filepath',physio_outdir);
        set_fpn = fullfile(physio_outdir,set_fname);

        outdata.data{ppaths_idx} = ensemble_add_data_struct_row(...
            outdata.data{ppaths_idx},'subject_id',subid,'session',isess,...
            'ensemble_id',sess.ensemble_id,'run',irun,'path',set_fpn,...
            'preprocessed',1);
        
        if SIG_CHECK
          sigs = {EEGf.chanlocs(:).labels};
          nsig = length(sigs);
          figure();
          for i=1:nsig
            sig = sigs{i};
            subplot(nsig,1,i);
            plot(EEGf.data(i,:));
            ylabel(sig);
          end % for i=1:nsig
          if isfield(defs,'figs') && isfield(defs.figs,'write2file') && ...
                  defs.figs.write2file
              figfname = fullfile(physio_outdir,sprintf(...
                  '%s_sess%d_run%d_physioFilteredSigCheck.ps',subid,isess,irun));
              print(figfname,defs.figs.printargs{:});
          end
        end % if SIG_CHECK
      end % if ~isempty(strmatch('filter
      
    end % for irun
    
  end % for isess
end % for isub=
