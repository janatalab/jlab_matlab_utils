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
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'realign_epi';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
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
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
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
check_vars = {'sinfo','epidata','pathdata'};
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
outdata.vars = epidata.vars;
ocol = set_var_col_const(outdata.vars);
for iv=1:length(outdata.vars)
  if isnumeric(epidata.data{iv})
    outdata.data{iv} = [];
  else
    outdata.data{iv} = {};
  end
end

% get flags
try RESLICE_EPI = defs.realign.RESLICE_EPI; catch RESLICE_EPI = 0; end
try USE_SPM = defs.realign.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.realign.USE_FSL; catch USE_FSL = 0; end

if ~USE_FSL && ~USE_SPM
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

%
% START SUBJECT LOOP
%

% generate subject masks
[sm,us] = make_mask_mtx(epidata.data{ocol.subject_id});

for is=1:length(us)
  % get subject
  subid = us{is};
  smask = sm(:,1);
  
  % get sinfo
  lsinfoidx = strmatch(subid,{sinfo.id});
  lsinfo = sinfo(lsinfoidx);
  
  % mask out other data
  subdata = epidata;
  for iv=1:length(outdata.vars)
    subdata.data{iv} = subdata.data{iv}(smask);
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
  
  [sessm,usess] = make_mask_mtx(subdata.data{ocol.session});
  for isess=1:length(usess)
    session = usess(isess);
    sesmask = sessm(:,isess);
    sessdata = subdata;
    for iv=1:length(outdata.vars)
      sessdata.data{iv} = sessdata.data{iv}(sesmask);
    end
    
    [runm,urun] = make_mask_mtx(sessdata.data{ocol.run});
    nruns = length(urun);

    sessinfo = lsinfo.sessinfo(session);
    
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
      for iv=1:length(outdata.vars)
        rundata.data{iv} = rundata.data{iv}(rm);
      end

      flist = rundata.data{ocol.path};
  
      if USE_SPM
        jobs{njob}.spatial{1}.realign{realign_est_idx}.estimate.data{irun} = ...
          cellstr(flist);
        if RESLICE_EPI
	      nf = size(flist,1);
	      jobs{njob}.spatial{realign_idx}.realign{realign_wrt_idx...
            }.write.data(end+1:end+nf) = cellstr(flist);
        end
      end

      % Do motion correction with FSL. Need to use a single 4D file. FSL
      % aligns things to the middle volume by default
      % GET run_outdir, fsl_mc_stub, subid
      if USE_FSL && USE_4D_NIFTI
        infname = flist;
	    outfname = fullfile(run_outdir,sprintf('%s_run%d_%s.nii', ...
	      subid, irun, fsl_mc_stub));
	  
        % Format the FSL command string
        fsl_str = sprintf('mcflirt -in %s -out %s -plots -rmsrel -rmsabs -report', ...
	      infname, outfname);
	    msg = sprintf('%s\n', fsl_str);
	    status = unix(fsl_str);  % execute the command
	    if status
	      msg = sprintf('UNIX command failed!!\n\n')
	      r = update_report(r,msg);
	      continue
        end
      end % if USE_FSL &

      % add flist contents to outdata
      for ifl=1:length(flist)
          outdata.data{ocol.subject_id} = ...
            [outdata.data{ocol.subject_id}; subid];
          outdata.data{ocol.session} = ...
            [outdata.data{ocol.session}; isess];
          outdata.data{ocol.ensemble_id} = ...
            [outdata.data{ocol.ensemble_id}; sessinfo.ensemble_id];
          outdata.data{ocol.run} = ...
            [outdata.data{ocol.run}; irun];
          outdata.data{ocol.path} = ...
            [outdata.data{ocol.path}; flist(ifl,:)];
      end
    end % for irun    
  end % for isess
end % for is

% Submit the SPM job stack
if RUN_SPM & ~isempty(jobs)
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
