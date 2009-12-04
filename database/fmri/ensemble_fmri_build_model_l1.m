function outdata = ensemble_fmri_build_model_l1(indata,defs)

% builds level 1 statistical model, optionally estimates model as well
% 
%   outdata = ensemble_fmri_build_model_l1(indata,defs)
% 
% This function builds and then estimates first-level general linear model
% analyses of fMRI data, using either SPM or FSL.
% 
% NOTE: this function expects that presentation logfiles exist for all runs
% that you want to analyze. If this data is not either provided through 
% 'pres_paths' in the indata, or findable on the disk using behav_indir,
% the function will not build the model.
% 
% REQUIRES
%   indata
%       sinfo
%       {'epi','realign_epi'}
%       paths
% 
% OPTIONAL INPUTS
%   indata
%       coplanar
%       modelspec - contains a pre-built model, if you want to build in one
%           step and then estimate in another step
%       physio_paths
%       pres_paths - if not found, will use 'behav_indir' from paths and
%           fmri_sess_paths, and will look for a typically-constructed
%           presentation data .mat file at that location. If not found, and
%           no pres_paths are given, the function will not run.
%   defs.fmri
%   defs.model
%   defs.model.permute_params
% 
% OUTPUT
%   sinfo
%   modelspec
%   resid
%   permutespec
% 
% FIXME: no permutation testing handled in FSL right now ....
% 
% FB 2008.08.27

global r

outdata = ensemble_init_data_struct();
outdata.type = 'build_model_l1';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
      switch indata{idata}.type
        case 'sinfo'
          sinfo = indata{idata}.data;
          proc_subs = {sinfo(:).id};
        case {'epi','realign_epi','resid'}
          epidata = indata{idata};
          epicol = set_var_col_const(epidata.vars);
        case 'paths'
          pathdata = indata{idata};
          pacol = set_var_col_const(pathdata.vars);
        case 'coplanar'
          coplanar = indata{idata};
          cocol = set_var_col_const(coplanar.vars);
        case 'modelspec'
          modelspec = indata{idata};
          mocol = set_var_col_const(modelspec.vars);
        case 'physio_paths'
          physio_paths = indata{idata};
          physcol = set_var_col_const(physio_paths.vars);
        case {'cardiac_epochs','pulse_epochs'}
          card_epochs = indata{idata};
          cepcol = set_var_col_const(card_epochs.vars);
        case {'gsr_epochs','scr_epochs'}
          gsr_epochs = indata{idata};
          gepcol = set_var_col_const(gsr_epochs.vars);
        case 'pres_paths'
          pres_paths = indata{idata};
          prescol = set_var_col_const(pres_paths.vars);
      end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'sinfo','epidata','pathdata'};
check_required_vars;

nsub_proc = length(proc_subs);

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
        sfilt.include.all.path_type = {'anal_outdir'};
        spathdata = ensemble_filter(lpathdata,sfilt);
        if length(spathdata.data{1}) == 1
          % one epi outdir, save outdata = epi_outdir
          outdata = spathdata.data{pacol.path}{1};
        else
          sfilt = pfilt;
          sfilt.include.all.path_type = {'sess_outdir'};
          spathdata = ensemble_filter(lpathdata,sfilt);
          if length(spathdata.data{1}) == 1;
            outdata = spathdata.data{pacol.path}{1};
          else
            sfilt = pfilt;
            sfilt.include.all.path_type = {'sub_outdir'};
            spathdata = ensemble_filter(lpathdata,sfilt);
            if length(spathdata.data{1}) == 1;
              outdata = spathdata.data{pacol.path}{1};            
            end
          end
        end
      end
    end
  end
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% get model
try curr_model = defs.model;
    model_id = curr_model.model_id;
catch
    msg = sprintf('Couldn''t load model from defs.model, aborting\n');
    r = update_report(r,msg);
    return
end

% get flags
try BUILD_MODEL = defs.build_model_l1.BUILD_MODEL; catch BUILD_MODEL=0; end % necessary??
try ESTIMATE_MODEL = defs.build_model_l1.ESTIMATE_MODEL; catch ESTIMATE_MODEL=0; end
try PLOT_DESMAT_COV = defs.build_model_l1.PLOT_DESMAT_COV; catch PLOT_DESMAT_COV=0; end
try USE_SPM = defs.build_model_l1.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.build_model_l1.USE_FSL; catch USE_FSL = 0; end
try SLICE_BASED = curr_model.slice_based; catch SLICE_BASED = 0; end

%
% get other vars
try modelperm_dir = params.build_model_l1.modelperm_dir; catch modelperm_dir = ''; end
try intensity_cutoff = params.build_model_l1.intensity_cutoff; catch intensity_cutoff = 10; end

% check flags/vars
if ~BUILD_MODEL && ~exist('modelspec','var')
    msg = sprintf('BUILD_MODEL not selected, but no model was specified');
    r = update_report(r,msg);
    return
end

if SLICE_BASED && (~USE_FSL || USE_SPM)
    msg = sprintf(['FSL is currently the only package supporting slice-'...
        'based modeling.']);
    r = update_report(r,msg);
    return
elseif SLICE_BASED
    msg = sprintf(['Slice-based modeling is not yet supported in this '...
        'script.']);
    r = update_report(r,msg);
    return
end

if isfield(curr_model,'combine_runs') && isnumeric(curr_model.combine_runs)...
    && ~isempty(curr_model.combine_runs) && curr_model.combine_runs 
  combine_runs = 1;
else
  combine_runs = 0;
end

if (~USE_FSL && ~USE_SPM) || (USE_FSL && USE_SPM)
  msg = sprintf(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  r = update_report(r,msg);
  return
end

if ESTIMATE_MODEL && ~BUILD_MODEL && USE_FSL && ~combine_runs
  error('this option is not yet supported')
end

if USE_SPM
  % Set stuff up for specifying an SPM job. Specify jobs on a subject level
  njob = 0;  % counter for number of jobs we are specifying
  jobs = {}; 

  % figure out if we're going to run the job
  RUN_SPM = defs.fmri.jobctl.run_spm;
elseif USE_FSL
  % set stuff up for specifying an FSL job.
  fsft = create_fsf; % fsf template
  if isfield(curr_model,'fsf_defs')
    fsft = set_fsf_defaults(fsft,curr_model.fsf_defs);
  else
    msg = sprintf('no fsf defaults provided');
    r = update_report(r,msg);
  end
end % if USE_SPM

if isfield(defs,'build_model_l1')
  if isfield(defs.build_model_l1,'USE_SPM_MOTPARAMS') && ...
          defs.build_model_l1.USE_SPM_MOTPARAMS
    USE_SPM_MOTPARAMS = 1;
    USE_FSL_MOTPARAMS = 0;
  elseif isfield(defs.build_model_l1,'USE_FSL_MOTPARAMS') && ...
          defs.build_model_l1.USE_FSL_MOTPARAMS
    USE_SPM_MOTPARAMS = 0;
    USE_FSL_MOTPARAMS = 1;
  else
    if USE_SPM
      USE_SPM_MOTPARAMS = 1;
      USE_FSL_MOTPARAMS = 0;
    elseif USE_FSL
      USE_SPM_MOTPARAMS = 0;
      USE_FSL_MOTPARAMS = 1;
    end
  end
end
exp_inroot = defs.paths.inroot;
exp_outroot = defs.paths.outroot;

% % % outdata

% sinfo
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% modelspec
outdata.vars = [outdata.vars 'modelspec'];
mod_idx = length(outdata.vars);
outdata.data{mod_idx} = ensemble_init_data_struct();
outdata.data{mod_idx}.type = 'modelspec';
outdata.data{mod_idx}.vars = {'subject_id','session','model_id','run','path'};
mcols = set_var_col_const(outdata.data{mod_idx}.vars);
outdata.data{mod_idx}.data{mcols.subject_id} = {};
outdata.data{mod_idx}.data{mcols.session} = [];
outdata.data{mod_idx}.data{mcols.model_id} = [];
outdata.data{mod_idx}.data{mcols.run} = [];
outdata.data{mod_idx}.data{mcols.path} = {};

% residuals
outdata.vars = [outdata.vars 'resid'];
res_idx = length(outdata.vars);
outdata.data{res_idx} = outdata.data{mod_idx};
outdata.data{res_idx}.type = 'resid';

if isfield(curr_model,'permute_params') && isstruct(curr_model.permute_params)
  % permutations
  outdata.vars = [outdata.vars 'permutespec'];
  perm_idx = length(outdata.vars);
  outdata.data{perm_idx} = outdata.data{mod_idx};
  outdata.data{perm_idx}.type = 'permutespec';
end

% 
% BEGIN SUBJECT LOOP
% 

for isub=1:nsub_proc
  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);
  
  % get paths
  spfilt.include.all.subject_id = {subid};
  spaths = ensemble_filter(pathdata,spfilt);
  
  didx = strmatch('sub_indir',spaths.data{pacol.path_type});
  if didx
    sub_indir = spaths.data{pacol.path}(didx);
  else
    sub_indir = fullfile(exp_inroot,subid);
    check_dir(sub_indir,1);
  end
  
  didx = strmatch('sub_outdir',spaths.data{pacol.path_type});
  if didx
    sub_outdir = spaths.data{pacol.path}(didx);
  else
    sub_outdir = fullfile(exp_outroot,subid);
    check_dir(sub_outdir,1);
  end
  
  % Here we need to set up output variables and file names that have to do
  % with subject specific things but which have to accommodate multiple sessions.
  if USE_SPM
    njob = njob+1;
    nspat = 0;  % counter  for number of spatial analyses within a job
    nstat = 0;  % counter for number of stats jobs
    ncon = 0;  % counter for number of contrast jobs
  end

  % Determine number of sessions for this subject
  nsess = length(sinfo(isub).sessinfo);
  
  %
  % START OF THE SESSION LOOP
  %
  
  for isess = 1:nsess
    sess = sinfo(isub).sessinfo(isess);
    
    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    else
      msg = sprintf('\n\nprocessing session %d\n\n',isess);
      r = update_report(r,msg);
    end
    
    % Deal with the path definitions
    fmri_sess_paths;
    
    if USE_SPM
      lanal_outdir = spm_outdir;
    elseif USE_FSL
      lanal_outdir = fsl_outdir;
    end
    
    % Figure out which version of the experiment applies
    exp_id = sess.exp_id;
    expidx = strmatch(exp_id,{defs.expinfo.id},'exact');
    if isempty(expidx)
      msg = sprintf('!!! Could not match experiment ID: %s\n', sess.exp_id);
      r = update_report(r,msg);
      continue
    end
    
    % Get the scanner protocol for this session
    protocol_idx = ...
	  strmatch(sess.protocol_id,{defs.fmri.protocol.id},'exact');
    protocol = defs.fmri.protocol(protocol_idx);
    
    if BUILD_MODEL

      % Load the Presentation files for this session
      if exist('pres_paths','var')
        presfilt.include.all.subject_id = {subid};
        presfilt.include.all.session = isess;
        lprespath = ensemble_filter(pres_paths,presfilt);
        pres_matfname = lprespath.data{prescol.path};
      end
      
      if ~exist('pres_matfname','var') || ~exist(pres_matfname,'file')
        pres_matfname = fullfile(behav_indir,sprintf('%s_sess%d_present.mat', subid, isess));
      end
      
      if exist(pres_matfname,'file')
        pres_info = load(pres_matfname);
        msg = sprintf('loaded presentation .mat file %s',pres_matfname);
        r = update_report(r,msg);
      else
        % no presentation data found, skip this session
        msg = sprintf(['presentation .mat file (%s) not found ... skipping '...
            '%s session %d\n'],pres_matfname,subid,isess);
        r = update_report(r,msg);
        continue
      end
      PL = set_pres_col_const(pres_info.vars);

      % Specify base model directory and make sure it exists
      model_outdir = fullfile(lanal_outdir, sprintf('model_%02d', model_id));
      check_dir(model_outdir);
      
      % Delete the contents of this directory
      unix_str = sprintf('rm -r %s', fullfile(model_outdir,'*'));
      status = unix(unix_str);

      if USE_SPM
        spmopts = defs.fmri.spm.opts(expidx).spmopts;
        build_model_idx_offset = nstat+1;
	
        % Initialize an fmri_spec using defaults
        fmri_spec = defs.fmri.spm.defaults.stats.fmri.fmri_spec;

        % Modify values accordingly
        fmri_spec.timing.RT = protocol.epi.tr;
	
        % Directory to output SPM.mat file to, Delete the contents of this directory
        unix_str = sprintf('rm %s', fullfile(model_outdir,'*'));
        status = unix(unix_str);

        fmri_spec.dir = {model_outdir};

        % Specify a mask file
        try mask_fname = curr_model.mask; 
          if isfield(curr_model,'maskmodel') && ~isempty(curr_model.maskmodel)
            mask_fname = fullfile(anal_outdir,'spm', ...
                sprintf('model_%02d', curr_model.maskmodel), mask_fname);
          end
        catch mask_fname = [];
        end
        
        if isempty(mask_fname)
          mask_fname = fullfile(anat_outdir,sprintf('sw%s_%s_mask.nii',...
              subid,sess.id));
        end
        if exist(mask_fname)
          fprintf('Using mask: %s\n', mask_fname);
          fmri_spec.mask = {mask_fname};
        end
          
        % Initialize a stats structure for this session
        nstat = nstat+1;

        jobs{njob}.stats{nstat}.fmri_spec = fmri_spec;
      end % if USE_SPM

      %
      % START OF RUN LOOP
      %

      % get # of runs
      if exist('epidata','var')
        % filter epidata by subject, session, run
        epiFilt = struct();
        epiFilt.include.all.subject_id = {subid};
        epiFilt.include.all.session = [isess];
        epiFilt.exclude.all.run = 0;
        sessdata = ensemble_filter(epidata,epiFilt);
        
        [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
        nruns = length(urun);
      else
        nruns = 0;
      end

      if USE_FSL && combine_runs
        % create cell to hold fsf structs
        fsf_structs = cell(nruns,1);
      end
      
      nvols = [];
      for irun = 1:nruns
        % get epi file list
        lmask = runm(:,irun);
        flist = sessdata.data{epicol.path}(lmask);
        
        % get presentation data for this run
        presfilt.include.all.RUN=irun;
        pinfo = ensemble_filter(pres_info,presfilt);
        
        if exist('physio_paths','var')
          lphys = ensemble_filter(physio_paths,presfilt);
          pinfo.physio_paths = lphys;
        end
        if exist('card_epochs','var')
          lcard = ensemble_filter(card_epochs,presfilt);
          pinfo.card_epochs = lcard;
        end
        if exist('gsr_epochs','var')
          lgsr = ensemble_filter(gsr_epochs,presfilt);
          pinfo.gsr_epochs = lgsr;
        end

        % Only execute analyses if this run exists or if we are dealing with a
        % model that is using residuals from a previous model, rather than the
        % EPI data
        try using_resid = strcmp(curr_model.srcdata,'residuals'); catch ...
	      using_resid = 0; end
      
        if ~isempty(flist) || using_resid

          % Check for congruency in number of pulses recorded by
          % Presentation and the number of volumes we actually have.

          % check # of pulses, etc
          pulsfilt.include.all.EVENT_TYPE = {'Pulse'};
          puls_info = ensemble_filter(pinfo,pulsfilt);
          pulse_times = puls_info.data{PL.RUN_REL_TIME};
          npulses = length(pulse_times);
          [cpp] = check_pulse_periods(pulse_times);

          predicted_nvol = npulses + cpp.data.nmiss;

          actual_nvol = length(flist);
          if actual_nvol == 1
            % is this a 4-d nifti file?
            unixstr = sprintf('fslval %s dim4',flist{1});
            [status,actual_nvol] = unix(unixstr);
            actual_nvol = str2num(actual_nvol);
          end
          nvols(irun) = actual_nvol;

          if actual_nvol < predicted_nvol
              msg = sprintf(['WARNING: Have fewer volumes (%d) than '...
                  'expected (%d)\n'], actual_nvol, predicted_nvol);
              r = update_report(r,msg);
          elseif actual_nvol > predicted_nvol
              msg = sprintf(['WARNING: Have more volumes (%d) than '...
                  'expected (%d)\n'], actual_nvol, predicted_nvol);
              r = update_report(r,msg);
          end

          estimated_tr = median(diff(pulse_times))/1000;
          if abs(diff([estimated_tr protocol.epi.tr])) > 0.02
            warning(sprintf(['WARNING: Stated TR (%f) and median pulse '...
                'period (%f), for %f/%f pulses, do not match!!!\n'],...
                protocol.epi.tr,estimated_tr,npulses,actual_nvol));
          end

          % Copy some scanning parameters to the run info structure
          epifstubfmt = defs.fmri.protocol(expidx).epi.fstubfmt;

% Add run specific stuff to the model structures
if USE_FSL
    
  %%%% FIXME: if multiple files in flist, concatenate
  epifname = flist{1};
  [fp,fn,fx] = fileparts(epifname);

  if ~isempty(regexp(fn,'\_bet','once'))
      fprintf(1,'brain extraction already applied, continuing\n');
  elseif exist(fullfile(fp,sprintf('%s_bet%s',fn,fx)),'file')
      fprintf(1,'brain extraction already applied, taking from disk\n');
      betfname = fullfile(fp,sprintf('%s_bet%s',fn,fx));
      nobetfname = epifname;
      epifname = betfname;
  else
      fprintf(1,'brain extraction: ');
      betfname = fullfile(fp,sprintf('%s_bet%s',fn,fx));
      betstr = sprintf('bet %s %s -F',epifname,betfname);
      status = unix(betstr);
      nobetfname = epifname;
      epifname = betfname;
      if status, error('error creating brain-extracted image'); end
      fprintf(1,'success!\n');
  end
  
  fsl_str = sprintf('fslval %s dim3', epifname);
  [status, nz] = unix(fsl_str);
  nslices = str2num(nz);

  fsl_str = sprintf('fslval %s dim4', epifname);
  [status, nvol] = unix(fsl_str);
  nvol = str2num(nvol);

  % Estimate the threshold intensity that we want to use for this run and
  % this subject
  fsl_str = sprintf('fslstats %s -M', epifname);
  [status, mean_intensity] = unix(fsl_str);
  mean_intensity = str2num(mean_intensity);
  thresh = mean_intensity*intensity_cutoff/100;

  % Calculate the mean of the EPI volume
  [fp,fn,fx] = fileparts(epifname);
  meanfname = fullfile(fp, sprintf('%s_mean%s',fn,fx));
  fsl_str = sprintf('fslmaths %s -Tmean %s', epifname, meanfname);
  status = unix(fsl_str);
  if status
    error('Failed to calculate mean of input EPI data file')
  end

  % create the fsf file
  fsf = fsft;
  fsf.npts = nvol;

  % generate EVs
  pinfo.evoutdir = fullfile(model_outdir,'evs'); check_dir(pinfo.evoutdir);
  pinfo.evstub = sprintf('%s_run%d_ev%%d.txt',subid,irun);
  pinfo.scanner.nslice_per_vol = nslices;
  pinfo.scanner.actual_nvol = nvol;
  pinfo.scanner.nvol = nvol;
  pinfo.scanner.dt = defs.fmri.spm.defaults.stats.fmri.fmri_spec.timing.fmri_t;
  pinfo.scanner.TR = protocol.epi.tr;
  pinfo.scanner.orig_nslices = protocol.epi.nslices;
  pinfo.protocol = protocol;
  pinfo.irun = irun;
  pinf.USE_FSL = 1;
  %%%%%%% HACK - FIXME: motionparam_fname should really either be a
  %%%%%%% filestub from params or better yet in its own output dataset from
  %%%%%%% ens_fmri_realign_epi
  if USE_FSL_MOTPARAMS
    pinfo.motionparam_fname = fullfile(fp,...
        sprintf('%s_%d_run%dmc.par',subid,sess.ensemble_id,irun));
  elseif USE_SPM_MOTPARAMS
    pinfo.motionparam_fname = fullfile(fp,...
        sprintf('rp_%s',sprintf(epifstubfmt,subid,irun,1,'txt')));
  end
  
  % fsf.ev = fmri_fsl_generate_evs(pinfo,curr_model,sess);
  fsf.evn = fmri_generate_regress(pinfo,curr_model,sess);
  nev = length(fsf.ev);

  fsf.evs_orig = nev;
  fsf.evs_real = nev;
  
  [fsf.ev.ortho] = deal(zeros(1,nev));
  
  fsf.tr = protocol.epi.tr;
  fsf.feat_files{1} = epifname;

  % Write the fsf file
  if combine_runs
    fsf_structs{irun} = fsf;
  else
    ev_names = {fsf.ev.name};
    contlist  = curr_model.contlist;
    Fcontlist = curr_model.Fcontlist;

    fsf = setup_fsl_con(fsf,ev_names,contlist,Fcontlist);

    run_outdir = fullfile(model_outdir,sprintf('run%d',irun));
    check_dir(run_outdir);
    fsf.fsldir = run_outdir;
    fsf.outputdir = run_outdir;
    
    % Retain a copy of the fsf structure so that we can use it during
    % subsequent statistical evaluation
    fsf_fname = fullfile(run_outdir, 'fsf.mat');
    save(fsf_fname, 'fsf');

    write_fsf(fsf);

    % Switch to the directory that we'll be evaluating the model in
    cd(run_outdir);

    % Check model integrity
    fsl_str = 'feat_model design';

    status = unix(fsl_str);
    if status
      msg = 'checking of model failed, SKIPPING';
      r = update_report(r,msg);
      continue
    end
    % Remove the old stats directory
    if exist('stats','dir')
      unix_str = 'rm -rf stats';
      unix(unix_str);
    end
    
    outdata.data{mod_idx} = ensemble_add_data_struct_row(outdata.data{mod_idx},...
        'subject_id',subid,'session',isess,'model_id',model_id,'run',irun,...
        'path',fullfile(run_outdir,'design.fsf'));

    if ESTIMATE_MODEL
      % Run the model
      fsl_str = sprintf('film_gls -rn ./stats -noest %s design.mat %1.4f',...
          epifname,thresh);
      status = unix(fsl_str);
      if status
        msg = 'FEAT run failed or was aborted';
        r = update_report(r,msg);
        continue
      end
    
      % Add the mean back into the residuals
      resid_fname =  fullfile(run_outdir,'stats','res4d.nii.gz');
      fsl_str = sprintf('fslmaths %s -add %s %s',...
          resid_fname,meanfname,resid_fname);
      status = unix(fsl_str);
      if status
        error('Failed to add mean back into the residuals')
      end

      outdata.data{res_idx} = ensemble_add_data_struct_row(outdata.data{res_idx},...
          'subject_id',subid,'session',isess,'model_id',model_id,'run',irun,...
          'path',resid_fname);
    end
    
  end % if combine_runs
      
elseif USE_SPM
	  
  % Add the list of input files
  spmsess.scans = cellstr(flist);
  
  pinfo.scanner.dt = defs.fmri.spm.defaults.stats.fmri.fmri_spec.timing.fmri_t;
  pinfo.scanner.actual_nvol = actual_nvol;
  pinfo.scanner.TR = protocol.epi.tr;
  pinfo.scanner.orig_nslices = protocol.epi.nslices;
  pinfo.irun=irun;
  pinfo.resp_mapping = sess.resp_mapping;
  pinfo.USE_SPM = 1;
  
  ppath = fileparts(flist{1});
  if USE_SPM_MOTPARAMS
    pinfo.motionparam_fname = fullfile(ppath,...
        sprintf('rp_%s',sprintf(epifstubfmt,subid,irun,1,'txt')));
  elseif USE_FSL_MOTPARAMS
    pinfo.motionparam_fname = fullfile(ppath,...
        sprintf('%s_%d_run%dmc.par',subid,sess.ensemble_id,irun));
  end
  pinfo.presfname = sess.pres.logfiles{irun,1};

  % if there is physio data for this subject/run, and if there is a physio
  % data struct that has been provided as input for this job, then extract
  % and attach the physo data to pinfo
  if isfield(sess,'physio')
    if exist('physio_paths','var')
      pf.include.all.subject_id = {subid};
      pf.include.all.session = isess;
      pf.include.all.run = urun(irun);
      lphys = ensembl_filter(physio_paths,pf);
      if ~isempty(lphys.data{physcol.path}) && ...
              exist(lphys.data{physcol.path}{1},'file')
        pinfo.physiofname = lphys.data{physcol.path}{1};
      else
        warning('physio file for run %d not found\n',urun(irun));
      end
    end
  end
  
  stat_idx = build_model_idx_offset;

  % Generate arrays specifying the conditions, onsets, durations
  lsess = fmri_generate_regress(pinfo,curr_model,sess);
  spmsess.cond = lsess.cond;
  spmsess.regress = lsess.regress;
    
  spmsess.multi = {''}; % {condinfo_fname}

  % Specify additional regressors
  % spmsess.regress = fmri_spm_generate_regress(pinfo,curr_model,sess);
  spmsess.multi_reg = {''};
  spmsess.hpf = curr_model.hpf;  % inf=no high-pass filtering

  % Add the current specification to the stack
  if irun == 1        
      % Calculate end of run time for this run
      end_of_run_time = actual_nvol*protocol.epi.tr;

      % remove condition onsets that occur after end_of_run_time
      for icond = 1:length(spmsess.cond)
        run_mask = spmsess.cond(icond).onset < end_of_run_time;
        if sum(run_mask) < length(spmsess.cond(icond).onset)
          spmsess.cond(icond).onset = ...
              spmsess.cond(icond).onset(run_mask);
          spmsess.cond(icond).duration = ...
              spmsess.cond(icond).duration(run_mask);
          if isfield(spmsess.cond(icond),'pmod') && ...
                  ~isempty(spmsess.cond(icond).pmod)
            % remove parametric modulators for condition onsets that occur
            % after end_of_run_time
            for imod=1:length(spmsess.cond(icond).pmod)
              spmsess.cond(icond).pmod(imod).param = ...
                spmsess.cond(icond).pmod(imod).param(run_mask);
            end
          end % if isfield(
        end % if sum(run_mask
      end % for icond
      jobs{njob}.stats{stat_idx}.fmri_spec.sess = spmsess;
  elseif ~combine_runs
      jobs{njob}.stats{stat_idx}.fmri_spec.sess(end+1) = spmsess;
  else
      % Make a local copy of the current session stack
      cs = jobs{njob}.stats{stat_idx}.fmri_spec.sess;
      init_num_scans = length(cs.scans);

      % Calculate offset (in sec) for this run
      offset = ...
        init_num_scans*jobs{njob}.stats{stat_idx}.fmri_spec.timing.RT;

      % Calculate end of run time for this run
      end_of_run_time = offset + (actual_nvol*protocol.epi.tr);
        
      % Concatenate list of scans
      cs.scans = [cs.scans; spmsess.scans];

      % Create a union of the condition lists
      allconds = union({cs.cond.name},{spmsess.cond.name});

      % Concatenate timing of conditions and copy condition information
      % in the event that a new condition has to be added to the list
      ncond = length(allconds);
      for icond = 1:ncond
		% Find the current condition in the original list
		curr_cond_idx = strmatch(allconds{icond},{cs.cond.name},'exact');
		new_cond_idx = strmatch(allconds{icond},{spmsess.cond.name},'exact');
        if ~isempty(new_cond_idx)
          % remove condition onsets that occur after end_of_run_time
          run_mask = (spmsess.cond(new_cond_idx).onset+offset) < end_of_run_time;
          if sum(run_mask) < length(spmsess.cond(new_cond_idx).onset)
            spmsess.cond(new_cond_idx).onset = ...
                spmsess.cond(new_cond_idx).onset(run_mask);
            spmsess.cond(new_cond_idx).duration = ...
                spmsess.cond(new_cond_idx).duration(run_mask);
            if isfield(spmsess.cond(new_cond_idx),'pmod') && ...
                    ~isempty(spmsess.cond(new_cond_idx).pmod)
              % remove parametric modulators for condition onsets that occur
              % after end_of_run_time
              for imod=1:length(spmsess.cond(new_cond_idx).pmod)
                spmsess.cond(new_cond_idx).pmod(imod).param = ...
                  spmsess.cond(new_cond_idx).pmod(imod).param(run_mask);
              end
            end
          end
          % concatenate conditions
          if ~isempty(curr_cond_idx)
            % concatenate condition onsets
		    cs.cond(curr_cond_idx).onset = ...
		        [cs.cond(curr_cond_idx).onset; ...
			  spmsess.cond(new_cond_idx).onset+offset];
            cs.cond(curr_cond_idx).duration = ...
                [cs.cond(curr_cond_idx).duration; ...
              spmsess.cond(new_cond_idx).duration];
            nnewmod = length(spmsess.cond(new_cond_idx).pmod);
            noldmod = length(cs.cond(curr_cond_idx).pmod);
            if nnewmod && noldmod
              % modulations previously added
              % concatenate parametric modulations
              for imod = 1:nnewmod
                oldmodidx = strmatch(...
                    spmsess.cond(new_cond_idx).pmod(imod).name,...
                    {cs.cond(curr_cond_idx).pmod(:).name});
                if oldmodidx
                  % this mod existed previously
                  % concatenate existing modulations
                  cs.cond(curr_cond_idx).pmod(oldmodidx).param = ...
                      [cs.cond(curr_cond_idx).pmod(oldmodidx).param; ...
                      spmsess.cond(new_cond_idx).pmod(imod).param];
                else
                  % this mod didn't exist previously
                  % add new modulation
                  noldmod = noldmod + 1;
                  cs.cond(curr_cond_idx).pmod(noldmod) = ...
                      spmsess.cond(new_cond_idx).pmod(imod);
                end % if oldmodidx
              end % for imod =
            elseif nnewmod
              % no previous modulations found
              % add new parametric modulation
              % warning, given that there were previously NO mods, now there
              % is at least One mod
              msg = sprintf(['new parametric modulations for previously '...
                  'unmodulated condition (%s), subject %s session %d '...
                  'run %d'],allconds{icond},isess,irun);
              r = update_report(r,msg);
              cs.cond(curr_cond_idx).pmod = spmsess.cond(new_cond_idx).pmod;
            end % if nnewmod && noldmod
          else
            % need to add the condition
  		    num_cond = length(cs.cond)+1;
  		    cs.cond(num_cond).name = spmsess.cond(new_cond_idx).name;
 		    cs.cond(num_cond).onset = ...
		        spmsess.cond(new_cond_idx).onset+offset;
		    cs.cond(num_cond).duration = ...
		        spmsess.cond(new_cond_idx).duration;
		    cs.cond(num_cond).tmod = spmsess.cond(new_cond_idx).tmod;
		    cs.cond(num_cond).pmod = spmsess.cond(new_cond_idx).pmod;
          end % if ~isempty(curr_cond_idx
        end % if ~isempty(new_cond_idx
	      end % for icond
	    
	      % Concatenate regressors
	      allregs = union({cs.regress.name},{spmsess.regress.name});
	      nreg = length(allregs);
	      for ireg = 1:nreg
		curr_name = allregs{ireg};
		% We need to handle some regressors as run specific, e.g. linear
		% trends, etc
		
		curr_reg_idx = strmatch(allregs{ireg},{cs.regress.name},'exact');
		new_reg_idx = ...
		    strmatch(allregs{ireg},{spmsess.regress.name},'exact');
	      
		if ~isempty(curr_reg_idx) && ~isempty(new_reg_idx)
		    cs.regress(curr_reg_idx).val = [cs.regress(curr_reg_idx).val; ...
			  spmsess.regress(new_reg_idx).val];
		elseif ~isempty(new_reg_idx)
		    num_regress = length(cs.regress)+1;
		    cs.regress(num_regress).name = ...
			spmsess.regress(new_reg_idx).name;
   		    cs.regress(num_regress).val = [zeros(init_num_scans,1); ...
			  spmsess.regress(new_reg_idx).val];
        end
      end % for ireg

      % fill out the end of regressors that weren't treated in this run
      for icr = 1:length(cs.regress)
          ncr = length(cs.regress(icr).val);
          nsc = length(cs.scans);
          if ncr < nsc
              cs.regress(icr).val = [cs.regress(icr).val; ...
                  zeros(nsc-ncr,1)];
          end
      end

      % Assign the modified structure back into the session structure
      jobs{njob}.stats{stat_idx}.fmri_spec.sess = cs;

  end % if irun == 1

end % if USE_FSL / elseif USE_SPM

      end % if exist_epi
    end % for irun

    % Add a constant for this run as a regressor if we are combining
    % regressors across runs
    if combine_runs
      if USE_SPM && ((isfield(curr_model,'srcdata') && ...
              ~strcmp(curr_model.srcdata,'residuals')) || ...
              ~isfield(curr_model,'srcdata'))

        stat_idx = build_model_idx_offset;

        cs = jobs{njob}.stats{stat_idx}.fmri_spec.sess;
        cwrun = 0;
        for irun = 2:length(nvols)
          constant_reg = zeros(length(cs.scans),1);
          constant_reg(sum(nvols(1:irun))-nvols(irun)+1:sum(nvols(1:irun))) = 1;
          const_reg_idx = length(cs.regress)+1;
          cs.regress(const_reg_idx).name = sprintf('run%d_constant', irun);
          cs.regress(const_reg_idx).val = constant_reg;
        end

        % Fill out any short regressors with zeros
        for ireg = 1:length(cs.regress)
          if length(cs.regress(ireg).val) < sum(nvols)
            cs.regress(ireg).val(sum(nvols)) = 0;
          end
        end
      
       jobs{njob}.stats{stat_idx}.fmri_spec.sess = cs;

        % save residual file name
        resid_fname =  fullfile(jobs{njob}.stats{stat_idx}.fmri_spec.dir{1},...
            'ResMS.hdr');

        outdata.data{res_idx} = ensemble_add_data_struct_row(outdata.data{res_idx},...
            'subject_id',subid,'session',isess,'model_id',model_id,'run',irun,...
            'path',resid_fname);
        
      elseif USE_FSL

        % use first run's fsf struct as template
        % concatenate EVs across runs (concatenate EVs of same name)
        % since linear EVs for each run have unique run-specific names, you
        % needn't handle these specially, but you must model run-mean for
        % each run separately.

        fsf = fsf_structs{1};
        fsf.outputdir = model_outdir;
        fsf.fsldir = model_outdir;

        evnames = {fsf.ev.name};
        nev = length(fsf.ev);
        desmat = zeros(fsf.npts,nev);
        
        % get EVs for the first run
        for iev=1:nev
          evfname = fsf.ev(iev).fname;
          desmat(:,iev) = cell2mat(loadtxt(evfname));
        end

        % get EVs for additional runs        
        for irun = 2:length(fsf_structs)
          lfsf = fsf_structs{irun};
          fsf.feat_files = [fsf.feat_files lfsf.feat_files];

          % row mask for current run
          rowmask = fsf.npts+1:fsf.npts+lfsf.npts;

          % add rows to the design matrix
          desmat(rowmask,:) = 0;
          
          % get number of EVs for the current run
          nlev = length(lfsf.ev);
          for iev=1:nlev
            lev = lfsf.ev(iev);

            % find the current EV in the existing EV names
            evidx = strmatch(lev.name,evnames,'exact');
            if isempty(evidx)
              % new regressor, add a column
              nev = nev+1;
              desmat(:,nev) = 0;

              % add name to evnames
              evnames = [evnames lev.name];
              
              % add new regressor values for this run
              desmat(rowmask,end) = cell2mat(loadtxt(lev.fname));

    	      % add ev struct to fsf
              fsf.ev(nev) = lev;
            elseif length(evidx) > 1
              error('more than one match found in evnames for EV %s',...
                  lev.name);
            else
              % current regressor, add to the column
              desmat(rowmask,evidx) = cell2mat(loadtxt(lev.fname));
            end % if isempty(evidx
          end % for iev=1:

          % add run mean regressor
          rmrname = sprintf('mean_run%d',irun);
          if isempty(strmatch(rmrname,evnames))
            nev = nev+1;
            desmat(:,nev) = 0;
            desmat(rowmask,end) = 1;
            evnames = [evnames rmrname];

            fsf.ev(nev) = create_ev;
            fsf.ev(nev).name = rmrname;
            fsf.ev(nev).shape = 2;
            %fsf.ev(nev) %% FIXME
          end % if isempty(strmatch(sprintf('mean_run%d
          fsf.npts = fsf.npts + lfsf.npts;
        end % for irun=2:

        % keep track of number of EVs in fsf struct
        fsf.evs_orig = nev;
        fsf.evs_real = nev;
        
        % save new EVs out to file
        ev_outdir = fullfile(fsf.fsldir,'evs');
        evstub = sprintf('%s_sess%d_ev%%d.txt',subid,isess);
        for iev=1:nev
          fsf.ev(iev).fname = fullfile(ev_outdir,sprintf(evstub,iev));
          regdata = desmat(:,iev);
          save(fsf.ev(iev).fname,'regdata','-ascii');
          %% HACK!!! FIXME
          % fix orthogonalization vectors for this ev
          fsf.ev(iev).ortho(length(fsf.ev(iev).ortho)+1:nev) = 0;
        end
        
        % set up contrasts
        contlist  = curr_model.contlist;
        Fcontlist = curr_model.Fcontlist;

        fsf = setup_fsl_con(fsf,evnames,contlist,Fcontlist);

        % merge EPIs?
    	epi_outdir = fileparts(fileparts(fsf.feat_files{1}));
        epifname = fullfile(epi_outdir,sprintf('%s_sess%d_concat_all',...
            subid,isess));
        fslstr = sprintf('fslmerge -t %s %s',epifname,...
            cell2str(fsf.feat_files,' '));
        status = unix(fslstr);
        if status
          error('error merging epi run data');
        end
        fsf.feat_files = {epifname};

        % REMOVE MEAN?

        % Retain a copy of the fsf structure so that we can use it during
        % subsequent statistical evaluation
        fsf_fname = fullfile(fsf.fsldir, 'fsf.mat');
        save(fsf_fname, 'fsf');

        write_fsf(fsf);

        % Switch to the directory that we'll be evaluating the model in
        cd(fsf.fsldir);

        % Check model integrity
        fsl_str = 'feat_model design';

        status = unix(fsl_str);
        if status
          msg = 'checking of model failed, SKIPPING';
          r = update_report(r,msg);
          continue
        end
        
        % Remove the old stats directory
        if exist('stats','dir')
          unix_str = 'rm -rf stats';
          unix(unix_str);
        end

        outdata.data{mod_idx} = ensemble_add_data_struct_row(outdata.data{mod_idx},...
            'subject_id',subid,'session',isess,'model_id',model_id,'run',irun,...
            'path',fullfile(fsf.fsldir,'design.fsf'));

        if ESTIMATE_MODEL
          % Run the model
          fsl_str = sprintf('film_gls -rn ./stats -noest %s design.mat %1.4f',...
              epifname,thresh);
          status = unix(fsl_str);
          if status
            msg = 'FEAT run failed or was aborted';
            r = update_report(r,msg);
            continue
          end

          % Add the mean back into the residuals
          resid_fname =  fullfile(fsf.fsldir,'stats','res4d.nii.gz');
          fsl_str = sprintf('fslmaths %s -add %s %s',...
              resid_fname,meanfname,resid_fname);
          status = unix(fsl_str);
          if status
            error('Failed to add mean back into the residuals')
          end

          outdata.data{res_idx} = ensemble_add_data_struct_row(outdata.data{res_idx},...
              'subject_id',subid,'session',isess,'model_id',model_id,'run',irun,...
              'path',resid_fname);
        end

	if PLOT_DESMAT_COV
	  msg = sprintf('Loading model info in: %s\n', fsf.fsldir);
      	  r = update_report(r,msg);
      
	  opts.printfig = 1;
      	  opts.title = sprintf('%s, session %d, model %02d', subid, isess, model_id);
      	  opts.figstub = fullfile(defs.paths.figpath, sprintf('descorr_model_%02d', model_id));
      	  if ~started_desmat_cov_file
            opts.append_str = '';
            started_desmat_cov_file = 1;
	  else
            opts.append_str = '-append';
          end
      	  figure(1), clf
      	  opts.figh = gcf;

      	  post_process_corrmat_fsl(fsf.fsldir,opts);
	end

      end % if USE_SPM && ((isfield(curr_model,'srcdata
    end % if combine_runs
    
    % Replace list of original file names with list of residual images
    if USE_SPM && isfield(curr_model,'srcdata') && ...
	  strcmp(curr_model.srcdata,'residuals')

      resid_flist = '';
      if isfield(curr_model,'srcdata') && ...
            strcmp(curr_model.srcdata,'residuals')
        srcdir = fullfile(lanal_outdir, sprintf('model_%02d', ...
	      curr_model.srcmodel),'resid');
        srcstub = sprintf('ResI*.img');
        resid_flist = get_spm_flist(srcdir,srcstub);
      end

      stat_idx = build_model_idx_offset;
      jobs{njob}.stats{stat_idx}.fmri_spec.sess.scans = cellstr(resid_flist);
    end

    if USE_SPM
      outmodpath = fullfile(lanal_outdir,sprintf('model_%02d',model_id),'SPM.mat');
    elseif USE_FSL
      outmodpath = fullfile(lanal_outdir,sprintf('model_%02d',model_id),'design.fsf');
    end

    if combine_runs
      outdata.data{mod_idx} = ensemble_add_data_struct_row(outdata.data{mod_idx},...
          'subject_id',subid,'session',isess,'model_id',model_id,...
          'path',outmodpath);
    end
    
end % if BUILD_MODEL
    
    if (ESTIMATE_MODEL || PLOT_DESMAT_COV) && USE_SPM
      if ~BUILD_MODEL
        mfilt.include.all.subject_id = {subid};
        mfilt.include.all.session = isess;
        mfilt.include.all.model_id = model_id;
        mdata = ensemble_filter(modelspec,mfilt);
        model_fname = mdata.data{mocol.path};
        model_dir = fileparts(model_fname);
      else
        model_dir = fullfile(lanal_outdir,sprintf('model_%02d', model_id));
        model_fname = fullfile(model_dir,'SPM.mat');
      end
    end

    if ESTIMATE_MODEL && USE_SPM
      nstat = nstat+1;
      jobs{njob}.stats{nstat}.fmri_est.spmmat = {model_fname};
      
      est_method = curr_model.est_method;
      switch est_method
        case {'Classical'}
          jobs{njob}.stats{nstat}.fmri_est.method.Classical = 1;
        case {'Bayesian'}
          jobs{njob}.stats{nstat}.fmri_est.method.Bayesian = curr_model.bayes;
        case {'Bayesian2'}
          jobs{njob}.stats{nstat}.fmri_est.method.Bayes2 = 1;
        otherwise
          msg = sprintf('ERROR: Unknown estimation method: %s\n', est_method);
          r = update_report(r,msg);
          continue
      end % switch est_method
    end % if ESTIMATE_MODEL & USE_SPM

    if PLOT_DESMAT_COV && USE_SPM
      msg = sprintf('Loading model info in: %s\n', model_fname);
      r = update_report(r,msg);
      
      opts.printfig = 1;
      opts.title = sprintf('%s, session %d, model %02d', subid, isess, model_id);
      opts.figstub = fullfile(defs.paths.figpath, sprintf('descorr_model_%02d', model_id));
      if ~started_desmat_cov_file
        opts.append_str = '';
        started_desmat_cov_file = 1;
      else
        opts.append_str = '-append';
      end
      figure(1), clf
      opts.figh = gcf;
      post_process_corrmat(model_fname, opts);
    end % if PLOT_DESMAT_COV && USE_SPM
  end % for isess
end % for isub=

% Submit the SPM job stack
if USE_SPM && RUN_SPM && ~isempty(jobs)
  % Save the job file so the we have a record of what we did
  tstamp = datenum(now);
  job_stub = sprintf('jobs_%s.mat', datestr(tstamp,30));
  check_dir(defs.paths.jobpath);
  job_fname = fullfile(defs.paths.jobpath, job_stub);
  msg = sprintf('Saving job info to: %s\n', job_fname);
  r = update_report(r,msg);
  save(job_fname,'jobs');
  
  warning off
  msg = sprintf('Launching SPM jobs ...\n');
  r = update_report(r,msg);
  spm_jobman('run',jobs);
  warning on
end  % if RUN_SPM

if exist('jobs','var') && ~isempty(jobs)
  r.data.jobs = jobs;
end

% 
% MODEL PERMUTATIONS
% 

if isfield(curr_model,'permute_params')
    
  % get permutation parameters
  pp = curr_model.permute_params;
  try maxiter=pp.maxiter; catch maxiter=150; end
  est_method = curr_model.est_method;
  switch est_method
    case {'Classical'}
      pjobs{1}.stats{1}.fmri_est.method.Classical = 1;
    case {'Bayesian'}
      pjobs{1}.stats{1}.fmri_est.method.Bayesian = curr_model.bayes;
    case {'Bayesian2'}
      pjobs{1}.stats{1}.fmri_est.method.Bayes2 = 1;
    otherwise
      msg = sprintf('ERROR: Unknown estimation method: %s\n', est_method);
      r = update_report(r,msg);
      return
  end % switch est_method
  stoprulefh = parse_fh(pp.stopping_rule_fun);
  
  % get SPM struct from previous job
  spmfname = fullfile(model_outdir,'SPM.mat');
  load(spmfname);
  origSPM = SPM;

  % calculate SPM
  permjob.spmmat = {spmfname};
  permjob.xCon = SPM.xCon;
  permjob.conspec = pp.conspec;
%   permjob.conspec.titlestr = pp.contrast.titlestr;
%   permjob.conspec.mask = []; % contrasts for masking
%   permjob.conspec.thresh = pp.contrast.thresh
%   permjob.conspec.extent = 1;
%   permjob.conspec.threshdesc = 'FWE';
%   permjob.conspec.contrasts = ...
%   strmatch('Tonreg_F',{SPM.xCon.name},'exact');
  
  [SPM,xSPM] = spm_getSPM(permjob);
  nvox = size(xSPM.XYZ,2);  
  pp.XYZ = xSPM.XYZ;
  pp.vol_dim = SPM.xVol.DIM';
  
  % load residuals from original job
  srcres_fname = fullfile(model_outdir,'ResMS.img');
  if ~exist(srcres_fname,'file'), error('residuals not found!'), end
  Vsrc = spm_vol(srcres_fname);  % map the filename
  fprintf('Reading source data: %s\n', srcres_fname);
  Ysrc = spm_read_vols(Vsrc); % read the data

  % create empty permute residual struct
  Vsrc = [];
  Yperm = zeros([size(Ysrc) 0]);
  
  % iterate until stopping rule or max iter
  done = false;
  iter = 0;
  permutations = []; % FIXME: save perms to make sure none are repeated ??
  while (done)
    iter = iter + 1;
    curr_dir = fullfile(model_outdir,sprintf('perm%04d',iter));
    check_dir(curr_dir);
    
    % permute, save data
    SPM = fmri_permute_desmat(origSPM,pp);
    outfname = fullfile(curr_dir,'SPM.mat');
    save(outfname,'SPM');
    
    % estimate model
    pjobs{1}.stats{1}.fmri_est.spmmat = {outfname};
    warning off
    msg = sprintf('Launching permutation %d ...\n',iter);
    r = update_report(r,msg);
    spm_jobman('run',pjobs);
    warning on
    
    % load new residuals
    permres_fname = fullfile(curr_dir,'ResMS.img');
    Vperm(iter) = spm_vol(permres_fname);
    Yperm(:,:,:,iter) = spm_read_vols(Vperm(iter));
    
    % estimate null distribution
    [Ymu,Ysig,Yprob,YpermZ] = fmri_permute_eval_dist(Ysrc,Yperm,xSPM.XYZ);
    
    % evaluate stopping rule
    pp.Ymu = Ymu;
    pp.Ysig = Ysig;
    pp.Yprob = Yprob;
    pp.YpermZ = YpermZ;
    [done,pp] = stoprulefh(pp);
    
    if iter == maxiter, done = true; end
  end

  % Write the probability volume to file
  Vprob = Vsrc;
  Vprob.fname = fullfile(srcpath,'PermProb.img');
  Vprob.descrip = 'Probability that ResMS is lower than permuted models';
  Vprob = rmfield(Vprob,'private');

  Vprob = spm_create_vol(Vprob);
  Vprob = spm_write_vol(Vprob,Yprob);  

  %%%% FIXME: make permutation report, including model stats and status
  

end

% % % % %
% % % % % Various sub-functions
% % % % %

function [con] = add_con_job(model_fname, continfo)
  global r
% We want to add a con job to the job stack.  The con structure
% requires two fields: .spmmat which specifies the SPM.mat file where the
% model and contrast information is stored and the consess structure
% which contains the actual contrast information.
con.spmmat = {model_fname};
consess = create_contrast_struct(model_fname,continfo);
con.consess = consess;
      
% Now we have to deal with cleaning up any existing contrast
% information in the model directory
      
% First we need to clean out the .Xcon structure in any existing
% SPM.mat file.
if exist(model_fname,'file')
  msg = sprintf('Clobbering xCon info in: %s\n', model_fname);
  r = update_report(r,msg);
  load(model_fname)
  SPM.xCon = {};
  save(model_fname,'SPM');
end
      
% Now remove any con, ess, spmF, and spmT files
cwd = pwd;
model_dir = fileparts(model_fname);
cd(model_dir);
files = {'^con_.{4}\..{3}$','^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};
for i=1:length(files)
  j = spm_select('List',pwd,files{i});
  for k=1:size(j,1)
    spm_unlink(deblank(j(k,:)));
  end
end
cd(cwd)

return
