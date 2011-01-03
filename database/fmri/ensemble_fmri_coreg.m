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

outdata = ensemble_init_data_struct();
outdata.type = 'coreg';

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
      case 'hires'
        hires = indata{idata};
        hicol = set_var_col_const(hires.vars);
      case 'coplanar'
        coplanar = indata{idata};
        cocol = set_var_col_const(coplanar.vars);
      case {'paths'}
        pathdata = indata{idata};
        pacol = set_var_col_const(pathdata.vars);
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

% get flags
try USE_SPM = defs.coreg.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.coreg.USE_FSL; catch USE_FSL = 0; end
try epi2coplanar_dof = defs.coreg.epi2coplanar_dof;
catch epi2coplanar_dof = 6; end % degrees of freedom

% get coreg list
coreg_list = defs.fmri.general.coreg_list;

% check for data for requested coregistrations
if ~isempty(strmatch('epi',{coreg_list{:,1:2}})) && ...
        ~exist('epidata','var')
  error('epi data requested in coregister list, but not provided');
end
if ~isempty(strmatch('coplanar',{coreg_list{:,1:2}})) && ...
        ~exist('coplanar','var')
  error('coplanars requested in coregister list, but not provided');
end
if ~isempty(strmatch('hires',{coreg_list{:,1:2}})) && ...
        ~exist('hires','var')
  error('hires data requested in coregister list, but not provided');
end

% outdata
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

if exist('epidata','var')
    % epidata
    outdata.vars = [outdata.vars 'epi'];
    epi_idx = length(outdata.vars);
    outdata.data{epi_idx} = ensemble_init_data_struct();
    outdata.data{epi_idx}.type='epi';
    outdata.data{epi_idx}.vars = epidata.vars;
    for iv=1:length(epidata.vars)
      if isnumeric(epidata.data{iv})
        outdata.data{epi_idx}.data{iv} = [];
      else
        outdata.data{epi_idx}.data{iv} = {};
      end
    end
end

if exist('coplanar','var')
    % coplanar data
    outdata.vars = [outdata.vars 'coplanar'];
    cop_idx = length(outdata.vars);
    outdata.data{cop_idx} = ensemble_init_data_struct();
    outdata.data{cop_idx}.type='coplanar';
    outdata.data{cop_idx}.vars = coplanar.vars;
    for iv=1:length(coplanar.vars)
      if isnumeric(coplanar.data{iv})
        outdata.data{cop_idx}.data{iv} = [];
      else
        outdata.data{cop_idx}.data{iv} = {};
      end
    end
end

if exist('hires','var')
    % coplanar data
    outdata.vars = [outdata.vars 'hires'];
    hires_idx = length(outdata.vars);
    outdata.data{hires_idx} = ensemble_init_data_struct();
    outdata.data{hires_idx}.type='hires';
    outdata.data{hires_idx}.vars = hires.vars;
    for iv=1:length(hires.vars)
      if isnumeric(hires.data{iv})
        outdata.data{hires_idx}.data{iv} = [];
      else
        outdata.data{hires_idx}.data{iv} = {};
      end
    end
end

if USE_FSL
    % transform matrix data
    outdata.vars = [outdata.vars 'xfm'];
    xfm_idx = length(outdata.vars);
    outdata.data{xfm_idx} = ensemble_init_data_struct();
    outdata.data{xfm_idx}.type='hires';
    outdata.data{xfm_idx}.vars = {'subject_id','session','ensemble_id',...
        'path_type','path'};
    xcol = set_var_col_const(outdata.data{xfm_idx}.vars);
    outdata.data{xfm_idx}.data{xcol.subject_id} = {};
    outdata.data{xfm_idx}.data{xcol.path} = {};
    outdata.data{xfm_idx}.data{xcol.path_type} = {};
    outdata.data{xfm_idx}.data{xcol.session} = [];
    outdata.data{xfm_idx}.data{xcol.ensemble_id} = [];
end

if (~USE_FSL && ~USE_SPM) || (USE_FSL && USE_SPM)
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

% set some initial paths
exp_inroot = defs.paths.inroot;
exp_outroot = defs.paths.outroot;

%
% START OF THE SUBJECT LOOP
%

for isub=1:nsub_proc

  % get subject
  subid = proc_subs{isub};
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

  lsinfo = sinfo(isub);
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
  end
  check_dir(sub_indir,1);
  
  didx = strmatch('sub_outdir',spaths.data{pacol.path_type});
  if didx
    sub_outdir = spaths.data{pacol.path}(didx);
  else
    sub_outdir = fullfile(exp_outroot,subid);
  end
  check_dir(sub_outdir,1);

  
  
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
  
  %
  % START OF THE SESSION LOOP
  %
  
  for isess = 1:nsess
    sessinfo = sess(isess);
    
    if ~sessinfo.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    end
    
    if isempty(sessinfo.use_epi_runs)
      msg = sprintf('\t\t\tno good runs, skipping session %d\n',isess);
      r = update_report(r,msg);
      continue
    end
    
    r = update_report(r,sprintf('\t\t\tsession %s\n', sessinfo.id));

    fmri_sess_paths;
    
    % Figure out which version of the experiment applies
    exp_id = sessinfo.exp_id;
    expidx = strmatch(exp_id,{defs.expinfo.id},'exact');
    if isempty(expidx)
      msg = sprintf('!!! Could not match experiment ID: %s\n', sessinfo.exp_id);
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
              r = update_report(r,msg);
              continue;
            end
            
            % get refvol from refrun
            epiFilt.include.all.subject_id = {subid};
            epiFilt.include.all.run = sessinfo.refrun;
            lepi = ensemble_filter(epidata,epiFilt);
            refvol = 1;

            srcfname = lepi.data{epicol.path}{refvol};
            if ~exist(srcfname,'file')
              msg = sprintf('WARNING: Desired source file does not exist: %s\n', srcfname);
              r = update_report(r,msg);
              continue
            end

          case 'coplanar'
            if ~exist('coplanar','var')
              msg = sprintf('WARNING: coplanar source data not found, %s sess %d',subid,isess);
              r = update_report(r,msg);
              continue;
            end
            srcfname = coplanar.data{cocol.path}{1};
            if isempty(strmatch(srcfname,outdata.data{cop_idx}.data{cocol.path}))
              outdata.data{cop_idx} = ensemble_add_data_struct_row(...
                  outdata.data{cop_idx},'subject_id',subid,'session',...
                  isess,'ensemble_id',sess.ensemble_id,'path',srcfname);
            end
          case 'hires'
            if ~exist('hires','var')
              msg = sprintf('WARNING: hires source data not found, %s sess %d',subid,isess);
              r = update_report(r,msg);
              continue;
            end
            srcfname = hires.data{hicol.path}{1};
            if isempty(strmatch(srcfname,outdata.data{hires_idx}.data{hicol.path}))
              outdata.data{hires_idx} = ensemble_add_data_struct_row(...
                  outdata.data{hires_idx},'subject_id',subid,'session',...
                  isess,'ensemble_id',sess.ensemble_id,'path',srcfname);
            end
          otherwise
            display(sprintf(['Do not know how to handle source type (%s).' ...
		      ' Skipping this coregistration\n'],src_type));
        end

        jobs{njob}.spatial{coreg_idx}.coreg{offset_idx}.estimate.source{1} = ...
	      srcfname;
	
        % Specify the reference (target) file
        switch ref_type
	      case 'epi'
            if ~exist('epidata','var')
              msg = sprintf('WARNING: epi target data not found, %s sess %d',subid,isess);
              r = update_report(r,msg);
              continue;
            end
	      case 'coplanar'
            if ~exist('coplanar','var')
              msg = sprintf('WARNING: coplanar target data not found, %s sess %d',subid,isess);
              r = update_report(r,msg);
              continue;
            end
            reffname = coplanar.data{cocol.path}{1};
          case 'hires'
            if ~exist('hires','var')
              msg = sprintf('WARNING: hires target data not found, %s sess %d',subid,isess);
              r = update_report(r,msg);
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
                r = update_report(r,msg);
                continue;
              end
              coplanar_fname = coplanar.data{cocol.path}{1};

              % Only add it if it exists
              if ~exist(coplanar_fname,'file')
                msg = sprintf('No coplanar image present: %s\n', coplanar_fname);
                r = update_report(r,msg);
              else
                other_flist{end+1} = coplanar_fname;
                if isempty(strmatch(coplanar_fname,...
                        outdata.data{cop_idx}.data{cocol.path}))
                  outdata.data{cop_idx} = ensemble_add_data_struct_row(...
                      outdata.data{cop_idx},'subject_id',subid,...
                      'session',isess,'ensemble_id',sess.ensemble_id,...
                      'path',coplanar_fname);
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
     
    elseif USE_FSL
      xfmdir = fullfile(sess_outdir, 'xfm');
      check_dir(xfmdir,0,1);
      refrunidx = sessinfo.use_epi_runs(1);
      filt = spfilt;
      filt.include.all.run = refrunidx;
      refrundata = ensemble_filter(epidata,filt);
      refrun_fname = refrundata.data{epicol.path}{1};
      volidx = get_middle_vol(refrun_fname);
      
      % Extract the timepoint
      pidxs = strfind(refrun_fname,'.');
      if isempty(pidxs)
        refvol_fname = sprintf('%s_t%d',refrun_fname,volidx);
      else
        refvol_fname = sprintf('%s_t%d%s',refrun_fname(1:pidxs(1)-1),...
            volidx,refrun_fname(pidxs(1):end));
      end
      fsl_str = sprintf('fslroi %s %s %d 1',refrun_fname,refvol_fname,volidx);
      display(fsl_str);
      unix(fsl_str);

      outdata.data{epi_idx} = ensemble_add_data_struct_row(...
          outdata.data{epi_idx},'subject_id',subid,'session',...
          isess,'ensemble_id',sess.ensemble_id,'run',refrunidx,...
          'path',refrun_fname);
          
    end % if USE_SPM
    
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
    end

    % iterate over runs
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
                outdata.data{epi_idx} = ensemble_add_data_struct_row(...
                    outdata.data{epi_idx},'subject_id',subid,'session',...
                    isess,'ensemble_id',sess.ensemble_id,'run',urun(irun),...
                    'path',flist(ifl,:));
              end
            end
          end
        elseif USE_FSL && ~isempty(strmatch('epi',{coreg_list{:,1}}))
          % if first run (reference run), skip
          if irun ==1, continue, end

          % extract middle volume
          infname = flist{1};
          targvol = get_middle_vol(infname);

          % Extract the timepoint for run 'irun'
          pidxs = strfind(infname,'.');
          if isempty(pidxs)
            targvol_fname = sprintf('%s_t%d',infname,targvol);
          else
            targvol_fname = sprintf('%s_t%d%s',infname(1:pidxs(1)-1),volidx,...
                infname(pidxs(1):end));
          end

          fsl_str = sprintf('fslroi %s %s %d 1',infname,targvol_fname,targvol);
          display(fsl_str);
          unix(fsl_str);

          % register the target volume to the reference volume
          xmat_fname = sprintf('%s_run%d_run%d.mat',subid,urun(irun),refrunidx);
          xmat_path = fullfile(xfmdir,xmat_fname);
          if exist(xmat_path,'file'), delete(xmat_path); end
          fsl_str = sprintf('flirt -in %s -ref %s -omat %s ',...
              targvol_fname,refvol_fname,xmat_path);
          display(fsl_str);
          unix(fsl_str);

          outdata.data{epi_idx} = ensemble_add_data_struct_row(...
              outdata.data{epi_idx},'subject_id',subid,'session',...
              isess,'ensemble_id',sess.ensemble_id,'run',urun(irun),...
              'path',infname);
          
          outdata.data{xfm_idx} = ensemble_add_data_struct_row(...
              outdata.data{xfm_idx},'subject_id',subid,'session',...
              isess,'ensemble_id',sess.ensemble_id,...
              'path_type',sprintf('run%d_2_run%d',urun(irun),refrunidx),...
              'path',xmat_path);
        end % if USE_SPM/elseif USE_FSL

      else
        msg = sprintf('run %d, empty flist, no residuals, SKIPPING',urun(irun));
        r = update_report(r,msg);
      end % if ~isempty(flist
    end % for irun
    
    if USE_FSL

      % set initial filter criteria
      filt.include.all.subject_id = {subid};
      filt.include.all.session = isess;
      
      % iterate over targets, coregister
      for ico=1:size(coreg_list,1)
        srcname = coreg_list{ico,1};
        refname = coreg_list{ico,2};
        
        if strcmp(srcname,refname)
          error(['source (%s) and reference (%s; coreg pair %d) are '...
              'the same!'],srcname,refname,ico);
        end

        % set source file name
        %%%% targvol_fname = srcfname
        switch srcname
          case 'epi'
            % get target epi fname
            srcfname = refvol_fname;
          case 'coplanar'
            % get coplanar_fname
            lsdata = ensemble_filter(coplanar,filt);
            srcfname = lsdata.data{cocol.path}{1};

            % save coplanar file to outdata
            outdata.data{cop_idx} = ensemble_add_data_struct_row(...
                outdata.data{cop_idx},'subject_id',subid,'session',...
                isess,'ensemble_id',sess.ensemble_id,'run',urun(irun),...
                'path',srcfname);
          case 'hires'
            % get hires_fname
            lsdata = ensemble_filter(hires,filt);
            srcfname = lsdata.data{hicol.path}{1};
            
            % save hires file to outdata
            outdata.data{hires_idx} = ensemble_add_data_struct_row(...
                outdata.data{hires_idx},'subject_id',subid,'session',...
                isess,'ensemble_id',sess.ensemble_id,'run',urun(irun),...
                'path',srcfname);
          otherwise
            error('unknown coregistration source: %s',srcname);
        end
        if ~exist(srcfname,'file') && ~exist([srcfname '.nii'],'file') ...
                && ~exist([srcfname '.nii.gz'],'file')
          error('source file %s not found',srcfname);
        end
        
        % set reference file name
        switch refname
          case 'epi'
            % if you get this error, check out the run loop starting in
            % this file, around lines 475-550
            error('epi coregistration handled differently');
          case 'coplanar'
            % get coplanar_fname
            lsdata = ensemble_filter(coplanar,filt);
            reffname = lsdata.data{cocol.path}{1};                
          case 'hires'
            % get hires_fname
            lsdata = ensemble_filter(hires,filt);
            reffname = lsdata.data{hicol.path}{1};    
          case 'template'
            % get template fname
            reffnamestub = sprintf('%s_template_path',srcname);
            reffname = defs.fmri.fsl.paths.(reffnamestub);
          otherwise
            error('unknown coregistration reference type: %s',refname);
        end
        if ~exist(reffname,'file') && ~exist([reffname '.nii'],'file') ...
                && ~exist([reffname '.nii.gz'],'file')
          error('reference file %s not found',reffname);
        end
        
        % Specify the name of the transformation matrix file
        xmat_fname = sprintf('%s_%s2%s.mat',subid,srcname,refname);
        xmat_path = fullfile(xfmdir,xmat_fname);
        
        % set degrees of freedom?
        if strcmp(refname,'coplanar') && strcmp(srcname,'epi')
          reg_dof = sprintf('-dof %d ',epi2coplanar_dof);
        else
          reg_dof = '';
        end
        
        % perform registration
        fsl_str = sprintf('flirt %s-in %s -ref %s -omat %s ',...
            reg_dof,srcfname,reffname,xmat_path);
        display(fsl_str);
        status = unix(fsl_str);
        if status
      	  error('ERROR: problem with coregistration of %s to %s',...
              srcname,refname);
        end
        
        % save xfm file to outdata
        outdata.data{xfm_idx} = ensemble_add_data_struct_row(...
            outdata.data{xfm_idx},'subject_id',subid,'session',...
            isess,'ensemble_id',sess.ensemble_id,...
            'path_type',sprintf('%s2%s',srcname,refname),...
            'path',xmat_path);
        
% % %       % Reslice hires
% % %       if RESLICE_HIRES
% % %     	unix_str = sprintf('rm %s', fullfile(anatdir,sprintf('r%s_hires.*', subid)));
% % % 	
% % %         reslice_fname = fullfile(anatdir,sprintf('r%s_hires.%s', subid, targ_fext_anat));
% % %         fsl_str = sprintf(['/usr/local/fsl/bin/flirt -in %s ' ...
% % %               '-ref %s -out %s -applyxfm -init %s'], hires_fname, template_fname, reslice_fname, xmat_fname);
% % %         fprintf(logfid,'%s\n', fsl_str);
% % %         status = unix(fsl_str);
% % %         if status
% % %           fprintf(logfid,['ERROR: problem with reslicing hires']);
% % %           error
% % %         end
% % %       end % if RESLICE_HIRES          
% % %           
% % %           %%%%%%%% NOTE: if epi to template ???
% % %           
% % %       xfm_fname = fullfile(xfmdir, sprintf('%s_run1_2_template.mat', subid));
% % %       hires2temp_xfm_fname = fullfile(xfmdir, sprintf('%s_hires2template.mat', subid));
% % %       epi2hires_xfm_fname = fullfile(xfmdir, sprintf('%s_run1_hires.mat', subid));
% % %       
% % %       fsl_str = sprintf(['/usr/local/fsl/bin/convert_xfm ' ...
% % % 	    '-omat %s -concat %s %s'], xfm_fname, hires2temp_xfm_fname, epi2hires_xfm_fname);
% % %       fprintf(logfid,'%s\n', fsl_str);
% % %       status = unix(fsl_str);
% % %       if status
% % % 	fprintf(logfid,['ERROR: problem with coregistration of run1 to' ...
% % % 	      ' template']);
% % % 	error
% % %       end
% % %       
% % %       
% % %       
% % %       

      end % for ico=1:size(coreg_list,1
    end % if USE_FSL
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
