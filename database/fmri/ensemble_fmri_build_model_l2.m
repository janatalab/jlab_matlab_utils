function outdata = ensemble_fmri_build_model_l2(indata,defs)

global r

outdata = ensemble_init_data_struct();
outdata.type = 'build_model_l2';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  switch indata(idata).type
    case 'modelspec'
      modelspec = indata(idata);
      mocol = set_var_col_const(modelspec.vars);
    case 'epi_mask_intersect'
      emidata = indata(idata);
      emicol = set_var_col_const(emidata.vars);
  end
end

if isfield(defs,'sinfo') && isstruct(defs.sinfo)
  sinfo = defs.sinfo;
  nsub_proc = length(sinfo(:));
end

% check for required vars, quit if they can't be found
check_vars = {'sinfo','modelspec','emidata'};
check_required_vars;

% outdata
% sinfo
outdata.vars = [outdata.vars; 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% modelspec
outdata.vars = [outdata.vars 'level2model'];
mod_idx = length(outdata.vars);
outdata.data{mod_idx} = ensemble_init_data_struct();
outdata.data{mod_idx}.type = 'level2model';
outdata.data{mod_idx}.vars = {'model_id','nsub','contrast','path'};
l2mocol = set_var_col_const(outdata.data{mod_idx}.vars);
outdata.data{mod_idx}.data{1} = [];
outdata.data{mod_idx}.data{2} = [];
outdata.data{mod_idx}.data{3} = {};
outdata.data{mod_idx}.data{4} = {};

% get model
try curr_model = defs.model;
    model_id = curr_model.model_id;
catch
    msg = sprintf('Couldn''t load model from defs.model, aborting\n');
    r = update_report(r,msg);
    return
end

try uses_permute = ~isempty(curr_model.permute.permute_type);
  uses_permute = 1;
catch
  uses_permute = 0;
end

if isfield(curr_model,'combine_runs') && isnumeric(curr_model.combine_runs)...
    && ~isempty(curr_model.combine_runs) && curr_model.combine_runs 
  combine_runs = 1;
else
  combine_runs = 0;
end

try BUILD_LEVEL2_MODEL = defs.build_model_l2.BUILD_LEVEL2_MODEL;
    catch BUILD_LEVEL2_MODEL = 0; end
try ESTIMATE_LEVEL2_MODEL = defs.build_model_l2.ESTIMATE_LEVEL2_MODEL;
    catch ESTIMATE_LEVEL2_MODEL = 0; end
try USE_SPM = defs.build_model_l2.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.build_model_l2.USE_FSL; catch USE_FSL = 0; end

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

exp_inroot = defs.paths.inroot;
exp_outroot = defs.paths.outroot;

% Deal with setting things up that operate across subjects, such as Level 2
% analyses
if BUILD_LEVEL2_MODEL | ESTIMATE_LEVEL2_MODEL
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
      % Specify where this analysis will be stored
      outpath = fullfile(defs.paths.analpath,'spm');
      check_dir(outpath);
      
      outpath = fullfile(outpath,'group');
      check_dir(outpath);
      
      outpath = fullfile(outpath, sprintf('model_%02d', model_id));
      check_dir(outpath);
      
      outpath = fullfile(outpath,level2(il).name);
      check_dir(outpath);
      
      model_fname = fullfile(outpath, 'SPM.mat');
      if exist(model_fname,'file')
        unixstr = sprintf('rm %s',model_fname);
        status = unix(unixstr);
        if status
          msg = sprintf(['error removing model file %s, job will fail if '...
              'we try to execute it in SPM, so quitting now\n'],model_fname);
          r = update_report(r,msg);
          return
        end
      end

      outdata.data{mod_idx}.data{l2mocol.model_id} = [...
          outdata.data{mod_idx}.data{l2mocol.model_id}; model_id];
      outdata.data{mod_idx}.data{l2mocol.nsub} = [...
          outdata.data{mod_idx}.data{l2mocol.nsub}; nsub_proc];
      outdata.data{mod_idx}.data{l2mocol.contrast} = [...
          outdata.data{mod_idx}.data{l2mocol.contrast}; level2(il).name];
      outdata.data{mod_idx}.data{l2mocol.path} = [...
          outdata.data{mod_idx}.data{l2mocol.path}; model_fname];
      
      % Load the structure template for the particular analysis we want to run
      if BUILD_LEVEL2_MODEL
        nstats=nstats+1;
        level2_idxs(il) = nstats;

        switch level2(il).type
          case 'one-sample-t'
            tmp = fd_t1;
        end

    	% Copy the structure template over
        jobs{njob}.stats{nstats} = tmp{1}.stats{1};
    	jobs{njob}.stats{nstats}.factorial_design.dir = {outpath};
	
        % Add the EPI intersection mask
        mask_fname = emidata.data{1};
        jobs{njob}.stats{nstats}.factorial_design.masking.em = mask_fname;
      end

      if ESTIMATE_LEVEL2_MODEL
        nstats=nstats+1;

        jobs{njob}.stats{nstats}.fmri_est.spmmat = {model_fname};
	
        if isfield(level2(il),'est_method')
            est_method = level2(il).est_method;
            switch est_method
              case {'Classical'}
                jobs{njob}.stats{nstats}.fmri_est.method.Classical = 1;
              case {'Bayesian'}
                jobs{njob}.stats{nstats}.fmri_est.method.Bayesian = curr_model.bayes;
              case {'Bayesian2'}
                jobs{njob}.stats{nstats}.fmri_est.method.Bayes2 = 1;
              otherwise
                msg = sprintf('ERROR: Unknown estimation method: %s\n', est_method);
                r = update_report(r,msg);
                continue
            end % switch est_method
        else
          jobs{njob}.stats{nstats}.fmri_est.method.Classical = 1;
        end
      end % if ESTIMATE_LEVEL2_MODEL
      
    end % for il = 1:nlevel2
  end % if ~level2
end % if BUILD_LEVEL2_MODEL

% 
% BEGIN SUBJECT LOOP
% x

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
    
    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    else
      msg = sprintf('\n\nprocessing session %d\n\n',isess);
      r = update_report(r,msg);
    end
        
    if BUILD_LEVEL2_MODEL & USE_SPM
      % get model spec
      modfilt.include.all.subject_id = {subid};
      modfilt.include.all.session = isess;
      lmod = ensemble_filter(modelspec,modfilt);
      
      if isempty(lmod.data{mocol.path}{1})
        msg = sprintf(['\nno model found for subject %s, session %d, '...
            'skipping!\n'],subid,isess);
        r = update_report(r,msg);
        continue
      end

      mpath = lmod.data{mocol.path}{1};
      model_indir = fileparts(mpath);
      
      % Get the xCon info from the SPM.mat file
      tmp = load(mpath);
      xCon = tmp.SPM.xCon;
      
      for il = 1:nlevel2
        cidx = level2_idxs(il);
        anal_type = ...
            fieldnames(jobs{level2_job_idx}.stats{cidx}.factorial_design.des);
        type_str = anal_type{1};
	
        switch type_str
          case 't1'
            scans = jobs{level2_job_idx}.stats{cidx}.factorial_design.des.(type_str).scans;
            conidx = strmatch(level2(il).src_contrast,{xCon.name});
            if ~isempty(conidx)
              % Recently patched versions of SPM5 seem to store the full path
              % and filename in the fname field, so check to see if this is the
              % case before building what would become an errorneous filename
              if exist(xCon(conidx).Vcon.fname,'file')
                scan_fname = xCon(conidx).Vcon.fname;
              else
                scan_fname = fullfile(model_indir, xCon(conidx).Vcon.fname);
              end
	      
              if iscell(scans)
                scans = [scans; scan_fname];
              else
                scans = {scan_fname};
              end
              jobs{level2_job_idx}.stats{cidx}.factorial_design.des.(type_str).scans = scans;
            end
        end % switch anal_type
      end % for il=1:nlevel2
    end % if BUILD_LEVEL2_MODEL & USE_SPM
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
