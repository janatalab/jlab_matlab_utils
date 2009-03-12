function outdata = ensemble_fmri_make4Dnifti(indata,defs)

% uses fslmerge to concatenate given images (.img, .nii) into a 4D nifti
%
% will merge both analyze and nifti (3d and 4d) into a new 4D nifti
% 
% REQUIRES
%   sinfo
%   epi data
%   path data
%   defs.output_dir_type = {epi_indir|epi_outdir|anat_indir|anat_outdir|...
%   	anal_indir|anal_outdir} ... essentially, any path type in path_data
%	default: epi_outdir
%   defs.CONCAT = {by_subject|by_run|all(default)}
% 
% 2009.03.05 FB

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'make4Dnifti';

r = init_results_struct;

r.type = 'make4Dnifti';  % Identify the type of this reporting instance
r.report_on_fly = 1;

% Parse out the input data
for idata = 1:length(indata)
  if isfield(indata{idata},'type')
    switch indata{idata}.type
      case 'sinfo'
        sinfo = indata{idata};
        sinfo = sinfo.data;
      case {'epi','realign_epi'}
        epidata = indata{idata};
        epicol = set_var_col_const(epidata.vars);
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
        outdata.data{epi_idx}.data{6} = {};
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
check_vars = {'sinfo','pathdata','epidata'};
check_required_vars;

% return an output directory for ensemble_jobman_parallel
% for this function, it will be the path_data type specified in
% defs.output_dir_type, defaulting to epi_outdir
if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
  if isfield(defs,'output_dir_type') && ischar(defs.output_dir_type) ...
        && ~isempty(defs.output_dir_type)
    odirtype = defs.output_dir_type;
  else
    odirtype = 'epi_outdir';
  end
  if exist('pathdata','var') && length(pathdata.data{1}) > 0
    if length(nsub_proc) == 1
      pfilt = struct();
      pfilt.include.all.subject_id = proc_subs;
      lpathdata = ensemble_filter(pathdata,pfilt);
      if ~isempty(lpathdata.data{1})
        sfilt = pfilt;
        sfilt.include.all.path_type = odirtype;
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
if ~isfield(defs,'CONCAT') || isempty(defs.CONCAT) ...
      || isempty(strmatch(defs.CONCAT,{'by_subject','by_run'}))
  CONCAT = 'all';
else
  CONCAT = defs.CONCAT;
end

switch CONCAT
    case 'all'
    otherwise



%
% START OF THE SUBJECT LOOP
%

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
    end
    
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
        
    end % for irun
  end % for isess
end % for isub=


end % switch CONCAT