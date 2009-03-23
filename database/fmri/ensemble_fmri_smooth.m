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
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
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
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
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
        sfilt.include.all.path_type = {'epi_outdir'};
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
  error('FSL not supported yet ...\n');
elseif ~USE_FSL && ~USE_SPM
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
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
      epiFilt.include.all.session = isess;
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
            outdata.data{epi_idx} = ensemble_add_data_struct_row(...
                outdata.data{epi_idx},'subject_id',subid,'session',isess,...
                'ensemble_id',sess.ensemble_id,'run',irun,'path',sfname);
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
if RUN_SPM && ~isempty(jobs)
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
