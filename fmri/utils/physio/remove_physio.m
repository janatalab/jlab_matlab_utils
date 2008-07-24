function remove_physio(proc,sinfo)
% remove_physio.m
%
% Removes variance associated with cardiac and respiration regressors on a
% slice-by-slice basis.
%
% Analysis steps to perform are specified in params which is a cell array of
% structures.
%
% proc{1}.type  {'eval_model','calc_stats','spat_reg'}
%
% Different analysis steps require additional control parameters which are
% specified in fields in the structure:
%  eval_model:
%    add_mean_to_resid - % Add the mean back to the estimated residuals. 
%    conlist
%    Fconlist

% 01/08/06 Petr Janata - apdated into function from remove_physio used for basic_states project

% Handle input arguments.  We need at least analysis parameters and subject
% info.
error(nargchk(2,2,nargin,'struct'))

% Read in the project globals
% bs_globals

nproc = length(proc);
defidx = [];
for iproc = 1:nproc
  if strcmp(proc{iproc}.type,'defs')
    defidx = iproc;
    rootdir = proc{iproc}.rootdir;
    verbose = proc{iproc}.verbose;
  end
end

stats_dir = 'stats';  % directory name within FEAT directory to store results in

logstub = '';
if ~isempty(logstub)
  logfile = fullfile(proc{defidx}.logdir, sprintf('%s_%s.txt',logstub,datestr(datenum(now),'yymmdd_HHMM')));
  logfid = fopen(logfile,'wt');
  fprintf('Writing logfile: %s\n', logfile);
else
  logfid = 1;
end

if proc{defidx}.use_mc moco_str = '_mc'; else moco_str = ''; end

start_dir = pwd;

nsub_proc = length(sinfo);
for isub = 1:nsub_proc
  subid = sinfo(isub).id;
  sdirs = check_subdirs(rootdir,subid,verbose);

  % Need to insert a session loop here
  
  % Load the physio information
  phys_fname = fullfile(sdirs.subdir, sprintf('%s_phys.mat', subid));
  fprintf(logfid,'Loading physio data from: %s\n', phys_fname);
  load(phys_fname)
  
  nruns = length(sinfo(isub).use_runs);
  for irun = 1:nruns
    run_id = sinfo(isub).use_runs(irun);
    
    % Make sure that run id in the phys structure matches the current run id
    if ~run_id == phys.ri(irun).id
      error(sprintf('Run ID mismatch: phys.ri.id=%d, run_id=%d', phys.ri(irun).id, run_id))
    end
    
    rdirs = check_rundirs(sdirs.physiodir, run_id,verbose);
    
    for iproc = 1:nproc
      cp = proc{iproc};
      switch cp.type
	case 'eval_model'
    
	  epifname = fullfile(sdirs.epidir, sprintf('Y_run%d%s_bet.nii.gz', run_id, moco_str));
	  % Check existence of the EPI file.  For some resting runs, data were not
	  % originally converted, so we want to skip over those
	  try
	    check_exist(epifname);
	  catch
	    fprintf(logfid,'WARNING: Wanted run %d but could not find %s\n', run_id, epifname);
	    continue
	  end
    
	  fsl_str = sprintf('avwval %s dim1', epifname);
	  [status, nx] = unix(fsl_str);
	  nx = str2num(nx);

	  fsl_str = sprintf('avwval %s dim2', epifname);
	  [status, ny] = unix(fsl_str);
	  ny = str2num(ny);

	  fsl_str = sprintf('avwval %s dim3', epifname);
	  [status, nz] = unix(fsl_str);
	  nslices = str2num(nz);
      
	  % Estimate the threshold intensity that we want to use for this run and
	  % this subject
	  fsl_str = sprintf('avwstats %s -M', epifname);
	  [status, mean_intensity] = unix(fsl_str);
	  mean_intensity = str2num(mean_intensity);
	  thresh = mean_intensity*cp.intensity_cutoff/100;

	  % Calculate the mean of the EPI volume
	  [fpath, fstub, fext] = fileparts(epifname);
	  meanfname = fullfile(fpath, sprintf('%s_mean.nii.gz', strtok(fstub,'.')));
	  fsl_str = sprintf('/usr/local/fsl/bin/avwmaths %s -Tmean %s', epifname, meanfname);
	  fprintf(logfid,'%s\n', fsl_str);
	  status = unix(fsl_str);
	  if status
	    error('Failed to calculate mean of input EPI data file')
	  end

	  for islice = 1:nslices
	    % Specify the output directory
	    featdir = fullfile(rdirs.rundir, sprintf('slice%02d.feat', islice));  % just the
	    % stub. 
	    check_dir(featdir);
	
	    funcfname = fullfile(sdirs.epidir, sprintf('slice%02d_data.nii.gz', islice));
	
	    % Create the fsf file
	    fsf = create_fsf;
	    fsf = set_fsf_defaults(fsf, 'remove_physio');
	
	    fsf.tr = cp.TR;
	    fsf.fsldir = featdir;  % determines where design.fsf will be written
	    fsf.outputdir = featdir;
	    fsf.feat_files{1} = funcfname;
	    fsf.npts = phys.ps.sinfo.nvol(run_id);
	
	    % Add the physio regressors
	    nsets = length(phys.ri(irun).ev_info);
	    nev = 0;
	    ev = create_ev;
	    ev_names = {};
	    
	    ev_info = phys.ri(irun).ev_info;
	    for iset = 1:nsets
	      nev_in_set = ev_info(iset).nev;
	      for iev = 1:nev_in_set
		cev = create_ev;
		cev.shape = 2;
		cev.fname = ev_info(iset).fnames{islice,iev};
	    
		nev = nev+1;
		ev(nev) = cev;
		ev_names{nev} = ev_info(iset).ev_names{iev};
	      end % for iev
	    end % for iset
	    fsf.ev = ev;
	
	    fsf.evs_orig = nev;
	    fsf.evs_real = nev;

	    %
	    % Deal with EV orthogonalization
	    % At the very least, we have to initialize these to zero
	    %
	    for iev = 1:nev
	      fsf.ev(iev).ortho = zeros(1,nev);
	    end
	
	    % Specify contrasts if we want to run the stats on the individual
	    % regressors.
	    % 
	    % NOTE: This has to be done at this stage rather than after packing all
	    % of the slice data into a single volume so that the proper template files
	    % are created for us to copy into the full volume analysis directory

	    fsf = setup_fsl_con(fsf, ev_names, cp.conlist, cp.Fconlist);

	    % Write the fsf file
	    write_fsf(fsf);
	
	    % Retain a copy of the fsf structure so that we can use it during
	    % subsequent statistical evaluation
	    fsf_fname = fullfile(featdir, 'fsf.mat');
	    save(fsf_fname, 'fsf');

	    % Switch to the directory that we'll be evaluating the model in
	    cd(featdir)
	
	    % Extract the relevant slice from the EPI data file
	    fsl_str = sprintf('avwroi %s %s 0 %d 0 %d %d 1 0 %d', epifname, funcfname, nx, ny, islice-1, phys.ps.sinfo.nvol(run_id));
	    fprintf(logfid,'%s\n', fsl_str);
	    status = unix(fsl_str);
	    if status
	      error('Failed to extract relevant slice from the EPI data file')
	    end

	    % Extract the relevant slice from the mean data file
	    slicemeanfname = fullfile(sdirs.epidir, sprintf('slice%02d_data_mean.nii.gz', islice));
	    fsl_str = sprintf('avwroi %s %s 0 %d 0 %d %d 1 0 1', meanfname, slicemeanfname, nx, ny, islice-1);
	    fprintf(logfid,'%s\n', fsl_str);
	    status = unix(fsl_str);
	    if status
	      error('Failed to extract relevant slice from the mean data file')
	    end
	
	    % Evaluate the model.
	    % Rather than calling FEAT which spawns a whole bunch of jobs, let's just
	    % run those components that we really want
	    
	    % Check model integrity
	    fsl_str = '/usr/local/fsl/bin/feat_model design';
	    
	    fprintf(logfid,'%s\n', fsl_str);
	    status = unix(fsl_str);
	    if status
	      error('Checking of model failed')
	    end
	
	    % Remove the old stats directory
	    if exist(stats_dir)
	      unix_str = sprintf('rm -rf %s', stats_dir);
	      unix(unix_str);
	    end
		
	    % Run the model
	    fsl_str = sprintf(['/usr/local/fsl/bin/film_gls -rn %s -noest' ...
		  ' %s design.mat %1.4f'], stats_dir, funcfname, thresh);
	    fprintf(logfid,'%s\n', fsl_str);
	    status = unix(fsl_str);
	    if status
	      error('FEAT run failed or was aborted')
	    end
	
	    if cp.add_mean_to_resid
	      % Add the mean back into the residuals
	      % NOTE: Adding the mean back has pros and cons.  It is good in the case that
	      % the residuals are passed into routines that extract the brain based on a
	      % proportion of the mean intensity in the image.  However, it is bad for some
	      % routines (like my cardiac-related BOLD estimation) that use the residuals
	      % image for masking.

	      resid_fname = fullfile(stats_dir, 'res4d.nii.gz');
	      fsl_str = sprintf('/usr/local/fsl/bin/avwmaths %s -add %s %s', resid_fname, slicemeanfname, resid_fname);
	      fprintf(logfid,'%s\n', fsl_str);
	      status = unix(fsl_str);
	      if status
		error('Failed to add mean back into the residuals')
	      end
	    end % if cp.add_mean_to_resid
	
	    % Remove the temporary functional data
	    delete(funcfname);
	    delete(slicemeanfname);
	
	    cd(start_dir)
	  end % for islice

	  % Pack all of the slice by slice residuals and parameter estimates into
	  % single volumes
            
	  % Get a list of the files that we want to pack together
	  stats_flist = dir(fullfile(rdirs.rundir,'slice01.feat/stats/*.nii.gz'));
	  nstats_images = length(stats_flist);
	  for iimg = 1:nstats_images
	    merge_list = {};
	    for islice = 1:nslices
	      merge_list{islice} = fullfile(rdirs.rundir, sprintf('slice%02d.feat', islice),'stats',stats_flist(iimg).name);
	    end
	
	    outfname = fullfile(rdirs.fullvolstatsdir,stats_flist(iimg).name);
	    fsl_str = sprintf('/usr/local/fsl/bin/avwmerge -z %s %s', outfname, cell2str(merge_list, ' '));
	    fprintf(logfid,'%s\n', fsl_str);
	    status = unix(fsl_str);
	    if status
	      error(sprintf('Failed to merge slice files into: %s', outfname))
	    end
	  end % for iimg = 1:nstats_images
	  % END OF EVAL_MODEL
    
	case 'calc_stats'  % Calculate significance maps using FEAT
	  fprintf(logfid, sprintf(['\nEvaluating statistical significance of the' ...
		' physiological variables in run %d\n'], run_id));
	  slice_featdir = fullfile(rdirs.rundir,'slice01.feat');

	  % Copy the mask that is associated with the original EPI volume
	  orig_maskfname = fullfile(sdirs.epidir, sprintf('Y_run%d%s_bet_mask.nii.gz', run_id, moco_str));
	  maskfname = fullfile(fullvoldir,'mask.nii.gz');
	  unix_str = sprintf('cp %s %s', orig_maskfname, maskfname);
	  status = unix(unix_str);
	  if status
	    error('Error copying mask file')
	  end
      
	  % Copy the degrees of freedom file from first slice analysis
	  unix_str = sprintf('cp %s %s', fullfile(slice_featdir,stats_dir,'dof'), rdirs.fullvolstatsdir);
	  status = unix(unix_str);
	  if status
	    error('Error copying dof file')
	  end

	  % Copy the original full volume time series
	  epifname = fullfile(sdirs.epidir, sprintf('Y_run%d%s_bet.nii.gz', run_id, moco_str));
	  funcfname = fullfile(rdirs.fullvoldir,'filtered_func_data.nii.gz');
	  unix_str = sprintf('cp %s %s', epifname, funcfname);
	  fprintf(logfid,'%s\n', unix_str);
	  status = unix(unix_str);
	  if status
	    error('Error copying original EPI data')
	  end
	
	  % Extract a single volume that we will use for overlays
	  example_funcfname = fullfile(rdirs.fullvoldir, 'example_func.nii.gz');
	  fsl_str = sprintf('/usr/local/fsl/bin/avwroi %s %s 0 1', funcfname, example_funcfname);
	  fprintf(logfid,'%s\n', fsl_str);
	  status = unix(fsl_str);
	  if status
	    error('Error extracting example functional volume')
	  end
	
	  % Load the fsf structure from slice 1 as a template
	  fsf_fname = fullfile(slice_featdir,'fsf.mat');
	  check_exist(fsf_fname);
	  load(fsf_fname)

	  % Modify the fsf structure accordingly
	  fsf = set_fsf_defaults(fsf, 'evaluate_physio');
	  fsf.outputdir = rdirs.fullvoldir;
	  fsf.fsldir = rdirs.fullvoldir;
	  fsf.feat_files{1} = rdirs.fullvoldir;  % FEAT directory
	
	  write_fsf(fsf);
	
	  % Save an fsf.mat file for posterity's sake
	  save(fullfile(rdirs.fullvoldir, 'fsf.mat'), 'fsf')

	  % Move to the 'FEAT' directory
	  cd(rdirs.fullvoldir)
	
	  fsl_str = sprintf('/usr/local/fsl/bin/feat %s', 'design.fsf');
	  fprintf(logfid,'%s\n', fsl_str);
	  status = unix(fsl_str);
	  if status
	    error('Problems running FEAT')
	  end
	  cd(start_dir)
	  
	  % end of calc_stats section
	  
	case 'spat_reg'  % perform spatial registration
	  % Load the existing fsf structure
	  fsf_fname = fullfile(rdirs.fullvoldir,'fsf.mat');
	  check_exist(fsf_fname);
	  load(fsf_fname)

	  fsf = set_fsf_defaults(fsf, 'coreg');
	  fsf.outputdir = rdirs.fullvoldir;
	  fsf.fsldir = rdirs.fullvoldir;
	  if ~iscell(fsf.feat_files)
	    fsf.feat_files = {rdirs.fullvoldir};
	  end
      
	  % Specify coplanar
	  coplanar_stub = sprintf('%s_coplanar', subid);
	  coplanar_fname = fullfile(sdirs.anatdir,sprintf('%s.nii.gz', coplanar_stub));
	  fsf.initial_highres_files{1} = coplanar_fname;
      
	  % Specify hires
	  hires_stub = sprintf('%s_hires', subid);
	  hires_fname = fullfile(sdirs.anatdir,sprintf('%s.nii.gz', hires_stub));
	  fsf.highres_files{1} = hires_fname;
      
	  write_fsf(fsf);

	  % Move to the 'FEAT' directory
	  cd(rdirs.fullvoldir)
      
	  fsl_str = sprintf('/usr/local/fsl/bin/feat %s', 'design.fsf');
	  fprintf(logfid,'%s\n', fsl_str);
	  status = unix(fsl_str);
	  if status
	    error('Problems running FEAT')
	  end
      
	  % Save an fsf.mat file for posterity's sake
	  save(fullfile(rdirs.fullvoldir, 'fsf.mat'), 'fsf')

	  cd(start_dir)
	  % end of spat_reg block
	otherwise
      end % switch proc.type
    end % for iproc=
  end % for irun=
end % for isub=
end % remove_physio

function sdirs = check_subdirs(rootdir,subid, verbose)
  subdir = fullfile(rootdir, subid);
  epidir = fullfile(subdir, 'epi');

  fsldir = fullfile(subdir, 'fsl');
  check_dir(fsldir, verbose);

  anatdir = fullfile(subdir, 'anat');
  check_dir(anatdir);
  
  physiodir = fullfile(subdir, 'physio_resid');
  check_dir(physiodir, verbose);

  % Populate the output struct
  sdirs = struct('subdir', subdir, 'epidir', epidir, 'fsldir', fsldir, ...
      'anatdir', anatdir, 'physiodir', physiodir);
end % check_subdirs

function rdirs = check_rundirs(rootdir, runid, verbose)
  rundir = fullfile(rootdir,sprintf('run%d', runid));
  check_dir(rundir, verbose);
  
  fullvoldir = fullfile(rundir, 'fullvol');
  check_dir(fullvoldir, verbose);
    
  fullvolstatsdir = fullfile(fullvoldir, 'stats');
  check_dir(fullvolstatsdir, verbose);
  
  % Populate the output struct
  rdirs = struct('rundir', rundir, 'fullvoldir', fullvoldir, ...
      'fullvolstatsdir', fullvolstatsdir);
end % check_rundirs
