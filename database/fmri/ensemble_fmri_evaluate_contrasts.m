function [outdata] = ensemble_fmri_evaluate_contrasts(indata,defs)

% evaluates contrasts for level1 models
% 
% FB 2008.08.27

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'evaluate contrasts';

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
        nsub_proc = length(proc_subs);
      case 'modelspec'
        modelspec = indata{idata};
        mocol = set_var_col_const(modelspec.vars);
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','modelspec','pathdata'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if exist('pathdata','var') && length(pathdata.data{1}) > 0
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
outdata.data{mod_idx}.vars = {'subject_id','session','model_id','path'};
outdata.data{mod_idx}.data{1} = {};
outdata.data{mod_idx}.data{2} = [];
outdata.data{mod_idx}.data{3} = [];
outdata.data{mod_idx}.data{4} = {};

% get model
try curr_model = defs.model;
    model_id = curr_model.model_id;
catch
    msg = sprintf('Couldn''t load model from defs.model, aborting\n');
    r = update_report(r,msg);
    return
end

% get flags
try USE_SPM = defs.build_model_l1.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.build_model_l1.USE_FSL; catch USE_FSL = 0; end

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
end % if USE_SPM

% 
% BEGIN SUBJECT LOOP
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
    else
      msg = sprintf('\n\nprocessing session %d\n\n',isess);
      r = update_report(r,msg);
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
      mfilt.include.all.subject_id = {subid};
      mfilt.include.all.session = isess;
      mfilt.include.all.model_id = model_id;
      mdata = ensemble_filter(modelspec,mfilt);
      model_fname = mdata.data{mocol.path}{1};
      
      spmopts = defs.fmri.spm.opts(expidx).spmopts;
      nstat = nstat+1;
      continfo = curr_model.continfo{1};
      jobs{njob}.stats{nstat}.con = add_con_job(model_fname, continfo);

      outdata.data{mod_idx}.data{1} = [outdata.data{mod_idx}.data{1}; subid];
      outdata.data{mod_idx}.data{2} = [outdata.data{mod_idx}.data{2}; isess];
      outdata.data{mod_idx}.data{3} = [outdata.data{mod_idx}.data{3}; model_id];
      outdata.data{mod_idx}.data{4} = [outdata.data{mod_idx}.data{4}; ...
        model_fname];
    end

  end % for isess
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


%
% Various sub-functions
%

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