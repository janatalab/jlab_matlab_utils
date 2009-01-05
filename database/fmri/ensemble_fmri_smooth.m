function outdata = ensemble_fmri_smooth(indata,defs)

% smooth data
% 
% REQUIRES
%   defs.smooth.USE_SPM
%   defs.smooth.USE_FSL
% 
% 2008.08.22 FB

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'smooth';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
      case {'epi','realign_epi'}
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
        outdata.vars = [outdata.vars 'epi'];
        epi_idx = length(outdata.vars);
        outdata.data{epi_idx} = ensemble_init_data_struct();
        outdata.data{epi_idx}.type='epi';
        outdata.data{epi_idx}.vars = epidata.vars;
        outdata.data{epi_idx}.data{1} = {};
        outdata.data{epi_idx}.data{2} = [];
        outdata.data{epi_idx}.data{3} = [];
        outdata.data{epi_idx}.data{4} = [];
        outdata.data{epi_idx}.data{5} = {};
    end
  end
end

if ~exist('sinfo','var')
  if isfield(defs,'sinfo')
    sinfo = defs.sinfo;
  end
end
proc_subs = {sinfo(:).id};
nsub_proc = length(proc_subs);

% check for required vars
check_vars = {'sinfo'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  outdata = defs.paths.outroot;
  if length(nsub_proc) == 1
    outdata = fullfile(outdata,proc_subs{1});
    if length(sinfo(1).sessinfo) == 1
      sdir = fullfile(outdata,'session1');
      if exist(sdir,'dir')
        outdata = sdir;
      end
    else
      % multiple sessions, save in the current outdir
    end
  else
    % multiple subjects, save in defs.paths.outroot
  end
  if ~exist(outdata,'dir'), outdata = ''; end
  return
end

% outdata
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% get flags
try USE_SPM = defs.smooth.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.smooth.USE_FSL; catch USE_FSL = 0; end

if USE_FSL && ~USE_SPM
  msg = sprintf('FSL not supported yet ...\n');
  r = update_report(r,msg);
  return
elseif ~USE_FSL && ~USE_SPM
  msg = sprintf(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  r = update_report(r,msg);
  return
end

% Set stuff up for specifying an SPM job.  Specify jobs on a subject level
if USE_SPM
  njob = 0;  % counter for number of jobs we are specifying
  jobs = {}; 

  % figure out if we're going to run the job
  RUN_SPM = defs.fmri.jobctl.run_spm;
end

%
% START OF THE SUBJECT LOOP
%

for isub=1:nsub_proc

  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

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
    end

    % Figure out which version of the experiment applies
    exp_id = sess.exp_id;
    expidx = strmatch(exp_id,{defs.expinfo.id},'exact');
    if isempty(expidx)
      msg = sprintf('!!! Could not match experiment ID: %s\n', sess.exp_id);
      r = update_report(r,msg);
      continue
    end
 
    if USE_SPM
      spmopts = defs.fmri.spm.opts(expidx).spmopts;
      nspat = nspat+1;
      smooth_idx = nspat;
      jobs{njob}.spatial{smooth_idx}.smooth.fwhm = spmopts.smooth_fwhm;
      smooth_flist = {};
      
      % SMOOTH COPLANAR?
    end

    %
    % START OF RUN LOOP
    %
    % get # of runs

    if exist('epidata','var')
      % filter epidata by subject, session, run
      epiFilt = struct();
      epiFilt.include.all.subject_id = {subid};
      epiFilt.include.all.session = [isess];
      sessdata = ensemble_filter(epidata,epiFilt);
        
      [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
      nruns = length(urun);
    else
      nruns = 0;
      warning('no epi images found ... nothing is being smoothed');
    end
      
    nvols = [];
    for irun = 1:nruns
      lmask = runm(:,irun);
      flist = sessdata.data{epicol.path}(lmask);

      % Only execute analyses if this run exists or if we are dealing with a
      % model that is using residuals from a previous model, rather than the
      % EPI data
      try using_resid = strcmp(curr_model.srcdata,'residuals'); catch ...
	    using_resid = 0; end
      
      if ~isempty(flist) || using_resid

	    if USE_SPM
          for ifl=1:length(flist)
            [fpath fname fext] = fileparts(flist{ifl});
            ofname = fullfile(fpath,sprintf('%s%s',fname,fext));
            smooth_flist = [smooth_flist; ofname];
            sfname = fullfile(fpath,sprintf('s%s%s',fname,fext));
            outdata.data{epi_idx}.data{epicol.subject_id} = ...
                [outdata.data{epi_idx}.data{epicol.subject_id}; subid];
            outdata.data{epi_idx}.data{epicol.session} = ...
                [outdata.data{epi_idx}.data{epicol.session}; isess];
            outdata.data{epi_idx}.data{epicol.ensemble_id} = ...
                [outdata.data{epi_idx}.data{epicol.ensemble_id}; sess.ensemble_id];
            outdata.data{epi_idx}.data{epicol.run} = ...
                [outdata.data{epi_idx}.data{epicol.run}; irun];
            outdata.data{epi_idx}.data{epicol.path} = ...
                [outdata.data{epi_idx}.data{epicol.path}; sfname];
          end

	      % Delete any previous smoothed images for good measure
          srcdir = fileparts(flist{1});
	      unix_str = sprintf('rm %s/sw%s*.??? >> /dev/null', srcdir, subid);
	      unix(unix_str);
          jobs{njob}.spatial{smooth_idx}.smooth.data = smooth_flist;
        end % if USE_SPM
      end % if exist_epi
    end % for irun=
  end % for isess=
end % for isub=

% Submit the SPM job stack
if RUN_SPM & ~isempty(jobs)
  % Save the job file so the we have a record of what we did
  tstamp = datenum(now);
  job_stub = sprintf('jobs_%s.mat', datestr(tstamp,30));
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
