function outdata = ensemble_fmri_ica(indata,defs)

% performs independent component analyses on the given volumes
%
% CURRENTLY only uses FSL's Melodic
% CURRENTLY expects 4D nifti volumes
% 
% REQUIRES
%   sinfo
%   epi data
%   path data
%   defs.ica.Z_SCORE - generates and uses z-scored images for analysis
%   defs.ica.CONCAT_BY_SUB - concatenates volumes for each subject
% 
% 2009.03.05 FB

global defaults r
VERBOSE = 1;

outdata = ensemble_init_data_struct();
outdata.type = 'ica';

r = init_results_struct;

r.type = 'fmri_ica';  % Identify the type of this reporting instance
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
        outdata.data{epi_idx}.type='epi';
        outdata.data{epi_idx}.vars = epidata.vars;
        outdata.data{epi_idx}.data{1} = {};
        outdata.data{epi_idx}.data{2} = [];
        outdata.data{epi_idx}.data{3} = [];
        outdata.data{epi_idx}.data{4} = [];
        outdata.data{epi_idx}.data{5} = {};
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
% for this function, it should be the analysis output directory
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
outdata.vars = [outdata.vars 'sinfo'];
sinfo_idx = length(outdata.vars);
outdata.data{sinfo_idx} = ensemble_init_data_struct();
outdata.data{sinfo_idx}.type = 'sinfo';
outdata.data{sinfo_idx}.data = sinfo;

% get flags
try USE_SPM = defs.ica.USE_SPM; catch USE_SPM = 0; end
try USE_FSL = defs.ica.USE_FSL; catch USE_FSL = 1; end
try CONCAT_BY_SUB = defs.ica.CONCAT_BY_SUB; catch CONCAT_BY_SUB = 0; end

if USE_SPM && ~USE_FSL
  error('SPM not supported yet ...\n');
elseif ~USE_FSL && ~USE_SPM
  error(['\t\tyou must specify either SPM or FSL to carry out '...
      'the analyses\n']);
end

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
