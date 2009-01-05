function outdata = ensemble_fmri_coreg(indata,defs)

% coregisters fmri data
% 
% REQUIRES
%   sinfo
%   hires data
%   coplanar data
%   epi data
% 
% FIXME: right now, it assumes it's only getting one subject, but it should
% really filter data by subject, session, and run where applicable to make
% sure it's respecting the subject information it's getting from sinfo
% 
% 2008.08.17 FB

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'coreg';

r = init_results_struct;

r.type = 'fmri_anal';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
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
  end
end

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
        sdir = fullfile(outdata,'epi');
        if exist(sdir,'dir')
          outdata = sdir;
        end
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

% get coreg list
coreg_list = defs.fmri.general.coreg_list;

% get flags
try USE_SPM = defs.realign.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.realign.USE_FSL; catch USE_FSL = 0; end

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

    nspat = nspat+1;
    num_coreg = size(coreg_list,1);
    jobs{njob}.spatial{nspat}.coreg = {};
    ncoreg = 0;  % counter for total number of cued coreg jobs
    coreg_idx = nspat;
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
      ncoreg = size(coreg_list,1);
      add_epi2other_flist = zeros(1,ncoreg);
      for ico = 1:ncoreg
        src_type = coreg_list{ico,1};
        ref_type = coreg_list{ico,2};

        % Calculate the offset into our coreg{} array
        offset_idx = (isub-1)*ncoreg+ico;

        % Specify the source file
        % NOTE: source file doesn't get returned in the outdata, since it
        % has not been transformed
        switch src_type
          case 'epi'
            if ~exist('epidata','var')
              msg = sprintf('WARNING: epi source data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
            
            % get refvol from refrun
            epiFilt.include.all.subject_id = {subid};
            epiFilt.include.all.run = sess.refrun;
            lepi = ensemble_filter(epidata,epiFilt);
            refvol = 1;

            srcfname = lepi.data{epicol.path}{refvol};
            if ~exist(srcfname)
              msg = sprintf('WARNING: Desired source file does not exist: %s\n', srcfname);
              r = update_report(r,msg);
              continue
            end

          case 'coplanar'
            if ~exist('coplanar','var')
              msg = sprintf('WARNING: coplanar source data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
            srcfname = coplanar.data{cocol.path}{1};
            if isempty(strmatch(srcfname,outdata.data{cop_idx}.data{cocol.path}))
              outdata.data{cop_idx}.data{cocol.subject_id} = ...
                [outdata.data{cop_idx}.data{cocol.subject_id}; subid];
              outdata.data{cop_idx}.data{cocol.session} = ...
                [outdata.data{cop_idx}.data{cocol.session}; isess];
              outdata.data{cop_idx}.data{cocol.ensemble_id} = ...
                [outdata.data{cop_idx}.data{cocol.ensemble_id}; sess.ensemble_id];
              outdata.data{cop_idx}.data{cocol.path} = ...
                [outdata.data{cop_idx}.data{cocol.path}; srcfname];
            end
          case 'hires'
            if ~exist('hires','var')
              msg = sprintf('WARNING: hires source data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
            srcfname = hires.data{hicol.path}{1};
            if isempty(strmatch(srcfname,outdata.data{hires_idx}.data{hicol.path}))
              outdata.data{hires_idx}.data{hicol.subject_id} = ...
                [outdata.data{hires_idx}.data{hicol.subject_id}; subid];
              outdata.data{hires_idx}.data{hicol.session} = ...
                [outdata.data{hires_idx}.data{hicol.session}; isess];
              outdata.data{hires_idx}.data{hicol.ensemble_id} = ...
                [outdata.data{hires_idx}.data{hicol.ensemble_id}; sess.ensemble_id];
              outdata.data{hires_idx}.data{hicol.path} = ...
                [outdata.data{hires_idx}.data{hicol.path}; srcfname];
            end
          otherwise
            fprintf(['Do not know how to handle source type (%s).' ...
		      ' Skipping this coregistration\n'], src_type);
        end

        jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.source{1} = ...
	      srcfname;
	
        % Specify the reference (target) file
        switch ref_type
	      case 'epi'
            if ~exist('epidata','var')
              msg = sprintf('WARNING: epi target data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
	      case 'coplanar'
            if ~exist('coplanar','var')
              msg = sprintf('WARNING: coplanar target data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
            reffname = coplanar.data{cocol.path}{1};
          case 'hires'
            if ~exist('hires','var')
              msg = sprintf('WARNING: hires target data not found, %s sess %d',subid,isess);
              r = update_report(r.msg);
              continue;
            end
            reffname = hires.data{hicol.path}{1};
          case {'template','canonicalT1'}
            reffname = fullfile(defs.fmri.spm.paths.canonical_dir,'avg152T1.nii');
          otherwise
            msg = sprintf(['Do not know how to handle source type (%s).' ...
              ' Skipping this coregistration\n'], src_type);
            r = update_report(r,msg);
        end
        jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.ref{1} = ...
	      reffname;
	
    	% Specify other files that we want to propagate the
        % transformations to
        other_flist = {};
        jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.other = {};
        other_types = coreg_list{ico,3};
        ntypes = length(other_types);
        for itype = 1:ntypes
          switch other_types{itype}
            case 'coplanar'
              if ~exist('coplanar','var')
                msg = sprintf('WARNING: coplanar propagation data not found, %s sess %d',subid,isess);
                r = update_report(r.msg);
                continue;
              end
              coplanar_fname = coplanar.data{cocol.path}{1};

              % Only add it if it exists
              if ~exist(coplanar_fname)
                msg = sprintf('No coplanar image present: %s\n', coplanar_fname);
                r = update_report(r,msg);
              else
                other_flist{end+1} = coplanar_fname;
                if isempty(strmatch(coplanar_fname,...
                        outdata.data{cop_idx}.data{cocol.path}))
                  outdata.data{cop_idx}.data{cocol.subject_id} = ...
                    [outdata.data{cop_idx}.data{cocol.subject_id}; subid];
                  outdata.data{cop_idx}.data{cocol.session} = ...
                    [outdata.data{cop_idx}.data{cocol.session}; isess];
                  outdata.data{cop_idx}.data{cocol.ensemble_id} = ...
                    [outdata.data{cop_idx}.data{cocol.ensemble_id}; sess.ensemble_id];
                  outdata.data{cop_idx}.data{cocol.path} = ...
                    [outdata.data{cop_idx}.data{cocol.path}; coplanar_fname];
                end
              end
            case 'epi'
              % Flag this coregistration job as one that should propage to
              % EPIs.  The file lists will be added in the run loop below.
              add_epi2other_flist(ico) = offset_idx;
            otherwise
              msg = sprintf('Unhandled other type: %s\n', other_types{itype});
              r = update_report(r,msg);
          end
	  	end % for itype
        
        jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.other = ...
	      other_flist;
	    jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.eoptions = ...
          spmopts.coreg_eopts;
      end % for ico=
     
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
      sessdata = ensemble_filter(epidata,epiFilt);
        
      [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
      nruns = length(urun);
    else
      nruns = 0;
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
      
      if ~isempty(flist) || using_resid;	

        if USE_SPM
          tmp_idxs = find(add_epi2other_flist);
          nf = size(flist,1);
          for ico = 1:length(tmp_idxs)
            offset_idx = add_epi2other_flist(tmp_idxs(ico));
            jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.other(end+1:end+nf) ...
              = cellstr(flist);
            for ifl=1:length(flist)
              if isempty(strmatch(flist(ifl,:),...
                      outdata.data{epi_idx}.data{epicol.path}))
                outdata.data{epi_idx}.data{epicol.subject_id} = ...
                  [outdata.data{epi_idx}.data{epicol.subject_id}; subid];
                outdata.data{epi_idx}.data{epicol.session} = ...
                  [outdata.data{epi_idx}.data{epicol.session}; isess];
                outdata.data{epi_idx}.data{epicol.ensemble_id} = ...
                  [outdata.data{epi_idx}.data{epicol.ensemble_id}; sess.ensemble_id];
                outdata.data{epi_idx}.data{epicol.run} = ...
                  [outdata.data{epi_idx}.data{epicol.run}; irun];
                outdata.data{epi_idx}.data{epicol.path} = ...
                  [outdata.data{epi_idx}.data{epicol.path}; flist(ifl,:)];
              end
            end
          end
        end
        
      end % if ~isempty(flist
    end % for irun
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
