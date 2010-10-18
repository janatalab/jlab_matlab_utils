function outdata = ensemble_fmri_realign_epi(indata,defs)

% realign EPI images
% 
% REQUIRES
%   defs.realign.RESLICE_EPI
%   defs.realign.USE_SPM
%   defs.realign.USE_FSL
% 
%   defs.fmri.jobctl.run_spm
%   defs.fmri.spm.opts(expidx).spmopts
% 
%   defs.expinfo.id
% 
% 2008-08-14 FB

global defaults r

outdata = ensemble_init_data_struct();
outdata.type = 'realign_epi';

r = init_results_struct;

r.type = 'fmri_realign';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'epi'
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case {'paths'}
        pathdata = indata{idata};
        pacol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','epidata','pathdata'};
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
        sfilt.include.all.path_type = {'epi_outdir'};
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

% outdata
outdata = ensemble_init_data_struct();
outdata.name = 'realign';
outdata.type = 'realign';

% epidata
outdata.vars = [outdata.vars 'epi'];
epi_idx = length(outdata.vars);
outdata.data{epi_idx} = ensemble_init_data_struct();
outdata.data{epi_idx}.type='epi';
outdata.data{epi_idx}.vars = epidata.vars;
for iv=1:length(outdata.data{epi_idx}.vars)
  if isnumeric(epidata.data{iv})
    outdata.data{epi_idx}.data{iv} = [];
  else
    outdata.data{epi_idx}.data{iv} = {};
  end
end

% sinfo
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% get flags
try RESLICE_EPI = defs.realign.RESLICE_EPI; catch RESLICE_EPI = 0; end
try USE_SPM = defs.realign.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.realign.USE_FSL; catch USE_FSL = 0; end

if (~USE_FSL && ~USE_SPM) || (USE_FSL && USE_SPM)
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  return
end

% Set stuff up for specifying an SPM job.  Specify jobs on a subject level
if USE_SPM
  % figure out if we're going to run the job
  RUN_SPM = defs.fmri.jobctl.run_spm;

  jobs = {}; 
  njob = 0;  % counter for number of jobs we are specifying
end

if USE_FSL
  fsl_mc_stub = defs.fmri.fsl.fsl_mc_stub;
end

exp_inroot = defs.paths.inroot;
exp_outroot = defs.paths.outroot;

%
% START SUBJECT LOOP
%

for is=1:nsub_proc

  % get subject
  subid = proc_subs{is};
  lsinfo = sinfo(is);
  sess = lsinfo.sessinfo;
  nsess = length(sess);

  if isempty([sess.use_epi_runs]) || ~any([sess.use_session])
    msg = sprintf('no good runs for subject %d (%s), SKIPPING',isub,subid);
    r = update_report(r,msg);
    continue;
  end
  
  % get paths
  spfilt.include.all.subject_id = {subid};
  spfilt.exclude.all.run = 0;
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
  
  % set up SPM job spec
  if USE_SPM
    njob = njob+1;

    nspat = 0;  % counter  for number of spatial analyses within a job
    nspat = nspat+1;
    % Set an index to keep track which of the spatial cell elements is for realign      
    realign_idx = nspat; 

    jobs{njob}.spatial{nspat}.realign = {};
    nrealign = 0;
  end

  for isess=1:nsess

    if ~sess(isess).use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    else
      msg = sprintf('\n\nprocessing session %d\n\n',isess);
      r = update_report(r,msg);
    end

    if isempty(sess(isess).use_epi_runs)
      msg = sprintf('\t\t\tno good runs, skipping session %d\n',isess);
      r = update_report(r,msg);
      continue
    end
    
    % set filt criteria
    filt = spfilt;
    filt.include.all.session = isess;
    
    sessdata = ensemble_filter(epidata,filt);
    
    fmri_sess_paths;

    [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
    nruns = length(urun);

    sessinfo = sess(isess);
    
    % Figure out which version of the experiment applies
    exp_id = sessinfo.exp_id;
    expidx = strmatch(exp_id,{defs.expinfo.id},'exact');
    if isempty(expidx)
      msg = sprintf('!!! Could not match experiment ID: %s\n', sess.exp_id);
      r = update_report(r,msg);
      continue
    end
        
    if USE_SPM
      spmopts = defs.fmri.spm.opts(expidx).spmopts;
      nrealign = nrealign+1;
      jobs{njob}.spatial{realign_idx}.realign{nrealign}.estimate.data = ...
        cell(1,nruns);
      jobs{njob}.spatial{realign_idx}.realign{nrealign}.estimate.eoptions = ...
        spmopts.realign_eopts;
      realign_est_idx = nrealign;
      if RESLICE_EPI
        nrealign = nrealign+1;
        jobs{njob}.spatial{realign_idx}.realign{nrealign}.write.roptions = ...
          spmopts.realign_wopts;
        jobs{njob}.spatial{realign_idx}.realign{nrealign}.write.data = {};
        realign_wrt_idx = nrealign;
      end
    end

    %
    % START OF RUN LOOP
    %

    for irun = 1:nruns
      run = urun(irun);
      rm  = runm(:,irun);
      
      rundata = sessdata;
      for iv=1:length(rundata.vars)
        rundata.data{iv} = rundata.data{iv}(rm);
      end

      flist = rundata.data{epicol.path};
  
      if USE_SPM
        jobs{njob}.spatial{1}.realign{realign_est_idx}.estimate.data{irun} = ...
          cellstr(flist);
        if RESLICE_EPI
	      nf = size(flist,1);
	      jobs{njob}.spatial{realign_idx}.realign{realign_wrt_idx...
            }.write.data(end+1:end+nf) = cellstr(flist);
        end
        
        % add flist contents to outdata
        for ifl=1:length(flist)
            outdata.data{epi_idx} = ensemble_add_data_struct_row(...
                outdata.data{epi_idx},'subject_id',subid,...
                'session',isess,'ensemble_id',sessinfo.ensemble_id,...
                'run',irun,'path',flist(ifl,:));
        end
      end

      % Do motion correction with FSL. Need to use a single 4D file. FSL
      % aligns things to the middle volume by default
      if USE_FSL
        if length(flist) > 1
          error('more than one EPI per run!? must use 4d nifti');
        end
        infname = flist{1};
        
        if ~exist(infname,'file') && ~exist([infname '.nii.gz'],'file')
          error('can not find nifti file %s',infname);
        end
        
        fslstr = sprintf('fslval %s dim4',infname);
        [status,nvols] = unix(fslstr);
        nvols = str2num(nvols);
        if nvols < 2
          error('one or less volumes in epi file %s',infname);
        end
        
        pidx = strfind(infname,'.');
        if isempty(pidx)
          outfname = sprintf('%s%s',infname,fsl_mc_stub);
        else
          outfname = sprintf('%s%s%s',infname(1:pidx(1)-1),fsl_mc_stub,...
              infname(pidx(1):end));
        end

        % Format the FSL command string
        fsl_str = sprintf(['mcflirt -in %s -out %s -plots -rmsrel '...
            '-rmsabs -report'],infname,outfname);
	    msg = sprintf('%s\n', fsl_str);
	    status = unix(fsl_str);  % execute the command
	    if status
	      msg = sprintf('UNIX command failed!!\n\n');
	      r = update_report(r,msg);
	      continue
        end
        
        % add flist contents to outdata
        outdata.data{epi_idx} = ensemble_add_data_struct_row(...
            outdata.data{epi_idx},'subject_id',subid,...
            'session',isess,'ensemble_id',sessinfo.ensemble_id,...
            'run',irun,'path',outfname);
      end % if USE_FSL &

    end % for irun    
  end % for isess
end % for is

% Submit the SPM job stack
if USE_SPM && RUN_SPM && ~isempty(jobs)
  % Save the job file so the we have a record of what we did
  tstamp = datenum(now);
  job_stub = sprintf('jobs_%s.mat', datestr(tstamp,30));
  job_fname = fullfile(defs.paths.jobpath, job_stub);
  check_dir(defs.paths.jobpath,1);
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
