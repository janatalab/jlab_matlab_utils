function outdata = ensemble_fmri_normalise(indata,defs)

% normalise images
% 
% REQUIRES
%   defs.normalise.RESAMP
%   defs.normalise.USE_SPM
%   defs.normalise.USE_FSL
% 
%   defs.fmri.jobctl.run_spm
%   defs.fmri.spm.opts(expidx).spmopts
%   defs.fmri.spm.paths.canonical_dir
%   defs.fmri.finfo.anat.hires_stub
%   defs.fmri.finfo.anat.coplanar_T1stub
% 
%   defs.expinfo.id
% 
% FIXME: should this script deal with placing output into a different
% directory than the input comes from?
% 
% 2008.08.22 FB

global defaults r
VERBOSE = 1;

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'normalise';

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
        case 'hires'
          hires = indata{idata};
          hicol = set_var_col_const(hires.vars);
          outdata.vars = [outdata.vars 'hires'];
          hires_idx = length(outdata.vars);
          outdata.data{hires_idx} = ensemble_init_data_struct();
          outdata.data{hires_idx}.type='hires';
          outdata.data{hires_idx}.vars = hires.vars;
          outdata.data{hires_idx}.data{1} = {};
          outdata.data{hires_idx}.data{2} = [];
          outdata.data{hires_idx}.data{3} = [];
          outdata.data{hires_idx}.data{4} = {};
        case 'coplanar'
          coplanar = indata{idata};
          cocol = set_var_col_const(coplanar.vars);
          outdata.vars = [outdata.vars 'coplanar'];
          cop_idx = length(outdata.vars);
          outdata.data{cop_idx} = ensemble_init_data_struct();
          outdata.data{cop_idx}.type='coplanar';
          outdata.data{cop_idx}.vars = coplanar.vars;
          outdata.data{cop_idx}.data{1} = {};
          outdata.data{cop_idx}.data{2} = [];
          outdata.data{cop_idx}.data{3} = [];
          outdata.data{cop_idx}.data{4} = {};
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
try RESAMP  = defs.normalise.RESAMP;  catch RESAMP = 0;  end
try USE_SPM = defs.normalise.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.normalise.USE_FSL; catch USE_FSL = 0; end

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

    nspat = nspat+1;
    jobs{njob}.spatial{nspat}.normalise = {};
    nnorm = 0;
    normalise_idx = nspat;
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
    
    r = update_report(r,sprintf('\t\t\t%s\n', sess.id));

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
      
      if exist('hires','var')
        % Create a job to estimate the parameters
        nnorm = nnorm+1;
        norm_est_idx = nnorm;

        % Add the source image
        hires_img = hires.data{hicol.path}{1};        
        anat_dir = fileparts(hires_img);

        jobs{njob}.spatial{normalise_idx}.normalise{norm_est_idx}.est.subj.source = ...
            {hires_img};
        jobs{njob}.spatial{normalise_idx}.normalise{norm_est_idx}.est.subj.wtsrc = ...
            {''};
        
        % Add the estimation options
        jobs{njob}.spatial{nspat}.normalise{norm_est_idx}.est.eoptions = ...
            spmopts.normalise_eopts;
      
        % Add the template image
        jobs{njob}.spatial{nspat}.normalise{norm_est_idx}.est.eoptions.template = ...
            {fullfile(defs.fmri.spm.paths.canonical_dir,'avg152T1.nii')};

        if RESAMP
          nnorm = nnorm+1;
      
          % Add the parameter file
          matfname =  fullfile(anat_dir, sprintf('%s_hires_sn.mat', subid));
          jobs{njob}.spatial{nspat}.normalise{nnorm}.write.subj.matname = {matfname};
      
          % Add the name of the original hires image
          jobs{njob}.spatial{nspat}.normalise{nnorm}.write.subj.resample = ...
            {hires_img};
      
          % Add the options for the hires image
          jobs{njob}.spatial{nspat}.normalise{nnorm}.write.roptions = ...
            spmopts.hires.normalise.ropts;
      
          % Check to see if the warped hires exists already and clobber it if
          % it does. Otherwise the job will die.
          hires_stub = defs.fmri.finfo.anat.hires_stub;
          tmpfname = fullfile(anat_dir, sprintf('w%s_%s', subid,hires_stub));
          [tmppath, tmpfstub] = fileparts(tmpfname);
          if exist(tmpfname)
            msg = sprintf('Deleting normalized hires ...\n');
            r = update_report(r,msg);
            unix_str = sprintf('rm %s/%s*', tmppath, tmpfstub);
            msg = sprintf('%s\n', unix_str);
            r = update_report(r,msg);
            unix(unix_str);
          end
          hires_outfname = tmpfname;
        else
          hires_outfname = hires_img;
        end % if RESAMP
        outdata.data{hires_idx} = ensemble_add_data_struct_row(...
            outdata.data{hires_idx},'subject_id',subid,'session',isess,...
            'ensemble_id',sess.ensemble_id,'path',hires_outfname);
      else
        msg = sprintf(['can''t find hires data for subject %s, quitting '...
            'normalise step'],subid);
        r = update_report(r,msg);
      end % if exist('hires

      norm_flist = {};  % Keep track of other images to propagate the
      % normalization info to

      if exist('coplanar','var')
          
        %%%%% FIXME: shouldn't the coplanar be entered here somewhere?
        nnorm = nnorm+1;
        % Copy the parameters from the previous jobs
        jobs{njob}.spatial{nspat}.normalise{nnorm} = ...
	    jobs{njob}.spatial{nspat}.normalise{nnorm-1};

        % Add the source image
        cop_img = coplanar.data{cocol.path}{1};        
        anat_dir = fileparts(hires_img);
        norm_flist = [norm_flist; cop_img];

        % Replace reslice options with the ones we want for EPI and coplanar
        jobs{njob}.spatial{nspat}.normalise{nnorm}.write.roptions = ...
	      spmopts.normalise_ropts;
      
        % Check to see if the warped coplanar exists already and clobber it if
        % it does. Otherwise the job will die.
        coplanar_T1_stub = defs.fmri.finfo.anat.coplanar_T1stub;
        tmpfname = fullfile(anat_dir, sprintf('w%s_%s', subid,coplanar_T1_stub));
        [tmppath, tmpfstub] = fileparts(tmpfname);
        if exist(tmpfname)
    	  msg = sprintf('Deleting normalized coplanar ...\n');
          r = update_report(r,msg);
          unix_str = sprintf('rm %s/%s*', tmppath, tmpfstub);
    	  msg = sprintf('%s\n', unix_str);
          r = update_report(r,msg);
    	  unix(unix_str);
        end
        outdata.data{cop_idx} = ensemble_add_data_struct_row(...
          outdata.data{cop_idx},'subject_id',subid,'session',isess,...
          'ensemble_id',sess.ensemble_id,'path',tmpfname);
      end % if exist('coplanar

    end % if USE_SPM

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
    end

    %
    % START OF RUN LOOP
    %
    
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
	
        %%%% DELETE OLD NORMALIZED FILES!!!!
        epi_dir = fileparts(flist{1});
        if ~isempty(get_spm_flist(epi_dir,sprintf('w%s*.img',subid)));
	      msg = sprintf('Deleting existing normalized files!!\n');
	      r = update_report(r,msg);
	      unix_str = sprintf('rm %s/w%s*', epi_dir,subid);
	      status = unix(unix_str);
        end
	
	    if USE_SPM
	      norm_flist = [norm_flist; cellstr(flist)];	  
        end

        for ifl=1:length(flist)
          [tmppath,tmpfname,tmpext] = fileparts(flist{ifl,:});
          outfname = fullfile(tmppath,sprintf('w%s%s',tmpfname,tmpext));
          outdata.data{epi_idx} = ensemble_add_data_struct_row(...
            outdata.data{epi_idx},'subject_id',subid,'session',isess,...
            'ensemble_id',sess.ensemble_id,'run',irun,'path',outfname);
        end
      end % if exist_epi
    end % for irun
    
    if USE_SPM
      jobs{njob}.spatial{normalise_idx}.normalise{nnorm}.write.subj.resample = ...
        norm_flist;      
    end    
  end % for isess
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
