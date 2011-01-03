function outdata = ensemble_fmri_split4Dnifti(indata,defs)

% uses fslsplit to split 4D images into 3D by subject, session, or run
%
% outdata = ensemble_fmri_split4Dnifti(indata,defs)
% 
% output locations:
%   a subject's epi_outdir - one subject, one session, multiple runs
%   a subject's epi_outdir, first session - one subject, multiple sessions
%   exproot/group_timeseries - multiple subjects
% 
% REQUIRES
%   sinfo
%   epi data
%   path data
%   defs.make4Dnifti.CONCAT = {by_subject|by_session|by_run|all(default)}
%   defs.make4Dnifti.bet_params = string of parameters to send to bet ...
%       default is '-m -f 0.4'
%   defs.make4Dnifti.tmpdir = path to directory to place temporary files
%   defs.make4Dnifti.name_stub = optional file name stub to be appended to
%       all files that are generated
% 
% RETURNS
%   4D nifti files in 'epi'
%   mean image of 4D nifti files in 'mean_epi'
%   std image of 4D nifti files in 'std_epi'
%   zscored 4D nifti files in 'zscore_epi'
%   bet masks for all 4D nifti files, in 'betmasks'
% 
% NOTE: this function will z-score and brain-extract all images that it
% receives, but will not perform any other pre-processing on the images.
% Therefore, you should motion-correct, realign, slice-time-correct,
% coregister, and perform any other pre-processing on the images before
% bringing them to this script. THIS SCRIPT ASSUMES ALL IMAGES HAVE BEEN
% COREGISTERED TO A COMMON SPACE.
% 
% NOTE: even if you request to concatenate all images across subjects,
% run-level .nii, session-level .nii, and subject-level .nii files will be
% generated.
%   run-level images will contain entries in all epidata columns.
%   session-level images will contain entries in the subject_id, session,
%       and path columns
%   subject-level images will contain entries in the subject_id and path
%       columns
%   across-subject images will only contain an entry in the path column.
%   where no 'entries' are found in columns, a [] or '' will be found,
%   depending on the class of the column contents.
% 
% NOTE: only supports by_run and by_session concatenation for now
% 
% FIXME: give option to take masks from other functions, instead of
% re-calculating masks within this function
% 
% FIXME: when generating masks, make sure to follow the same convention for
% output data structs as ensemble_fmri_mask_from_{mean|norm}_epi
% 
% NOTE: assumes you're making a 4D image file from EPIs, returns paths to
% 4d images in a data struct named 'epi'
% 
% 2009.03.05 FB

global r

outdata = ensemble_init_data_struct();
outdata.type = 'split4Dnifti';

r = init_results_struct;

r.type = 'split4Dnifti';  % Identify the type of this reporting instance
r.report_on_fly = 1;

dstr = datestr(now(),30);

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
        proc_subs = {sinfo(:).id};
        nsub_proc = length(proc_subs);
      case {'epi','realign_epi','epidata','resid'}
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
      case {'paths'}
        pathdata = indata{idata};
        pcol = set_var_col_const(pathdata.vars);
    end
  end
end

% check for required vars
check_vars = {'sinfo','pathdata','epidata'};
check_required_vars;

% return an output directory for ensemble_jobman_parallel
% for this function, it will be the path_data type specified in
% defs.output_dir_type, defaulting to epi_outdir
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  % get output directory
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
          % one output directory, save outdata = odirtype path
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
        end % if length(spatndata.data{1
      end % if ~isempty(lpathdata
    else
      if isfield(defs,'paths') && isfield(defs.paths,'grouptime') && ...
              ~isempty(defs.paths.grouptime) && exist(defs.paths.grouptime)
        outdata = defs.paths.grouptime;
      end
    end % if length(nsub_proc
  end % if exist('pathdata
  if ~exist('outdata','var') || ~exist(outdata,'dir'), outdata = ''; end
  return
end

% sinfo output data struct
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% epi output data struct
outdata.vars = [outdata.vars 'epi'];
epi_idx = length(outdata.vars);
outdata.data{epi_idx} = ensemble_init_data_struct();
outdata.data{epi_idx}.type=epidata.type;
outdata.data{epi_idx}.vars = epidata.vars;
outdata.data{epi_idx}.data{1} = {};
outdata.data{epi_idx}.data{2} = [];
outdata.data{epi_idx}.data{3} = [];
outdata.data{epi_idx}.data{4} = [];
outdata.data{epi_idx}.data{5} = {};

%
% START OF THE SUBJECT LOOP
%

sub_files = {};
sub_masks = {};

for isub=1:nsub_proc

  subid = sinfo(isub).id;
  msg = sprintf('\t\tPROCESSING SUBJECT (%d/%d): %s\n', isub, nsub_proc,subid);
  r = update_report(r,msg);

  % get subject epis
  sfilt.include.all.subject_id = {subid};
  sdata = ensemble_filter(epidata,sfilt);

  % get subject paths
  spdata = ensemble_filter(pathdata,sfilt);
  
  % Determine number of sessions for this subject
  nsess = length(sinfo(isub).sessinfo);

  %
  % START OF THE SESSION LOOP
  %

  sess_files = {};
  sess_masks = {};
  
  for isess = 1:nsess
    sess = sinfo(isub).sessinfo(isess);            

    if ~sess.use_session
      msg = sprintf('\t\t\tSkipping session %d\n', isess);
      r = update_report(r,msg);
      continue
    end

    % get session epis
    sessfilt.include.all.session = isess;
    sessdata = ensemble_filter(sdata,sessfilt);
    
    % get session paths
    sesspdata = ensemble_filter(spdata,sessfilt);

    %
    % START OF RUN LOOP
    %

    [runm,urun] = make_mask_mtx(sessdata.data{epicol.run});
    nruns = length(urun);

    run_files = {};
    run_masks = {};
    
    for irun = 1:nruns
      lrun = urun(irun);
      runfilt.include.all.run = lrun;
      rundata = ensemble_filter(sessdata,runfilt);

      flist = rundata.data{epicol.path};
      if length(flist) > 1
        error('more than one 4d nifti provided for sub %s, run %d',...
            subid,lrun);
      end
      
      runfilt.include.all.path_type = {'run_outdir'};
      pdata = ensemble_filter(pathdata,runfilt);
      outpath = pdata.data{pcol.path}{1};
      cwd = pwd;
      cd(outpath);
      
      out_basename = sprintf('split%s_epi_r%02d_',subid,irun);
      
      unixstr = sprintf('fslsplit %s %s -t',flist{1},out_basename);
      status = unix(unixstr);
      if status
        error('problem splitting %s',flist{1});
      end
      
      files = dir(fullfile(outpath,[out_basename '*']));

      cd(cwd);
      names = {files.name};
      for ii=1:length(names)
        lname = fullfile(outpath,names{ii});
        outdata.data{epi_idx} = ensemble_add_data_struct_row(...
            outdata.data{epi_idx},'subject_id',subid,'session',...
            isess,'ensemble_id',sess.ensemble_id,'run',lrun,'path',lname);
      end

    end % for irun=1:
  end % for isess=1:
end % for isub=1:
