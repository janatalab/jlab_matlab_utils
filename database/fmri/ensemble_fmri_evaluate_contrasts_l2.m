function outdata = ensemble_fmri_evaluate_contrasts_l2(indata,defs)

global r

outdata = ensemble_init_data_struct();
outdata.type = 'eval_contrasts_l2';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'level2model'
        modelspec = indata{idata};
        mocol = set_var_col_const(modelspec.vars);
    end
  end
end

% check for required vars, quit if they can't be found
check_vars = {'modelspec'};
check_required_vars;

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(defs,'paths') && isfield(defs.paths,'analpath')
    outdata = defs.paths.analpath;
    check_dir(outdata);
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

try USE_SPM = defs.evaluate_contrasts_l2.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.evaluate_contrasts_l2.USE_FSL; catch USE_FSL = 0; end

if USE_FSL && ~USE_SPM
  error('FSL not supported yet ...\n');
  return
elseif ~USE_FSL && ~USE_SPM
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
  return
end

% Set stuff up for specifying an SPM job.  Specify jobs on a subject level
if USE_SPM
  njob = 0;  % counter for number of jobs we are specifying
  jobs = {}; 

  % figure out if we're going to run the job
  RUN_SPM = defs.fmri.jobctl.run_spm;

  % Deal with setting things up that operate across subjects, such as Level 2
  % analyses
  level2 = curr_model.level2;
  nlevel2 = length(level2);

  if ~nlevel2
    msg = sprintf('\tNo Level 2 model specified for model: %s\n', curr_model.name);
    r = update_report(r,msg);
    BUILD_LEVEL2_MODEL = 0;
  else
    % Initialize a stats job for all of the level 2 analyses
    njob = njob+1;
    level2_job_idx = njob;
    nstats = 0;
    level2_idxs = [];
    
    for il = 1:nlevel2

      % find contrast
      modfilt.include.all.contrast = {level2(il).name};
      mdata = ensemble_filter(modelspec,modfilt);

      mfname = mdata.data{mocol.path}{1};
      
      if ~exist(mfname,'file')
        msg = sprintf('model file %d missing, skipping contrast %s\n\n',...
            mfname,level2(il).name);
        r = update_report(r,msg);
        continue
      else
        nstats=nstats+1;
        sprintf('adding con job\n')
        jobs{njob}.stats{nstats}.con = ...
            add_con_job(mfname,level2(il).continfo);
      end % if EVALUATE_LEVE2_CONTRASTS

    end % for il = 1:nlevel2
  end % if ~level2
end % if USE_SPM

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