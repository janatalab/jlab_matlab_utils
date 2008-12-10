% 
%  Specifications for batch processing data for simulate_batch
%
%  This can be viewed as a template batch specification file, though in the
%  simulate context it is being used only to create design matrices
%

%
%  GLOBAL variables
%

global subject_ids proc_subs use_model_proto sinfo model_action
global rootpath dataroot batch_root
global group_analyses

global SCAN_OFFSET                % has to be specified in get_subject_info.m
                                  % number of volumes at the start of
				  % each run that were "hidden" from
				  % view, i.e. excluded from the analysis
					
% Regressor globals
global REGRESS_MOTION_CORRECTION_PARAMS ADD_LINEAR_TREND_REGRESS ADD_RT_REGRESS
REGRESS_MOTION_CORRECTION_PARAMS = 0;
ADD_LINEAR_TREND_REGRESS = 0;
ADD_RT_REGRESS =0;

global COREG_WITH_MUTUAL_INFO		% Perform coregistration using mutual
                                        % information method rather than the
                                        % standard method.  This hopefully
                                        % overcomes some of the distortions due
                                        % to the FIL EPI template.

rootpath = dataroot; % dataroot is specified in simulate_start.m

%
%  STRUCTURE CONTAINING LIST OF ANALYSES TO PERFORM
%
%  The analyses structure is created dynamically based on the list of desired
%  analyses
%

SET_GENERAL_DEFS = 1;
SET_F_EQUAL1 = 1;

%
% PREPROCESSING OPTIONS
%

REALIGN_EPI = 0;
CREATE_MEAN_EPI = 0;

COREG_WITH_MUTUAL_INFO = 0;
COREG_MEAN_EPI = 0;
COREG_FIRST_EPI = 0;

COREG_COPLANAR_2_HIRES = 0;

COREG_HIRES_2_TEMPLATE = 0;
RESLICE_EPI = 0;

NORM_COPLANAR_2_TEMPLATE = 0;

NORM_EPI_2_TEMPLATE = 0;
NORM_EPI_2_TEMPLATE_CUST_BB = 0;
NORM_EPI_2_TEMPLATE_SPM99_DEF_BB = 0; %USE this
NORM_EPI_2_TEMPLATE_SPM99_BB4 = 0;

NORM_HIRES_2_TEMPLATE = 0;

SMOOTH_rEPI = 0;
SMOOTH_nEPI = 0;

%
% MODELLING OPTIONS
%

RUN_MODEL_INDIV = 0;			% Create and evaluate model of single
                                        % subject, unnormalized data
RUN_MODEL_NORM = 1;			% Create and evaluate model of single
                                        % subject, normalized data (for group analyses)
use_model_proto = 1;  
model_action = 'specify';

COMPUTE_CONTRASTS_INDIV = 0;
COMPUTE_CONTRASTS_NORM = 0;

COMPUTE_GROUP = 0;   % tw implemented _ only for t-contrasts
group_analyses = {'RC-UC_hrf', 'RD-UD_hrf', 'RC-RD_hrf', 'UC-UD_hrf'};

THRESH_SUM = 0; THRESHSUM_CONTRAST=2;

%
%  SUBJECT INFO
%
subject_ids = cellstr(char(sinfo(:).id))

nsub = length(subject_ids);
proc_subs = [1]			% Indices of subjects to process


%
%  LIST OF ANALYSIS TYPES
%

% This is a list of the analysis types supported by SPM99
type = {'defaults_edit','model','contrasts','headers',...
'means','realign','coreg','normalize','smooth','display','threshsum', 'stats'};

%
%  GROUP ANALYSIS PARAMETERS
%

anal_list = strvcat(type);

% SPM99 expects to find a structure called 'analyses' that contains the 4
% fields specified in the following line of code
analyses = struct('type',[],'index',[],'work_dir',[],'mfile',[]);

na = length(analyses.type);

%
%  LIST OF WORK DIRECTORIES
%

work_dir = {rootpath, fullfile(rootpath, 'group_deriv'), ...
      };

for isub = 1:nsub
  work_dir{end+1} = fullfile(rootpath, subject_ids{isub}, 'single');
end
 
for isub = 1:nsub
  work_dir{end+1} = fullfile(rootpath, subject_ids{isub}, 'stats');
end

%
%  SPECIFY ANALYSES TO RUN
%

if SET_GENERAL_DEFS
  disp('Will set general defaults ...')
  na = na+1;
  analyses.type(na) = strmatch('defaults_edit', anal_list,'exact');
  analyses.index(na) = 1;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 1;   
end

if REALIGN_EPI
  disp('Will realign EPI runs ...')
  na = na+1;
  analyses.type(na) = strmatch('realign', anal_list,'exact');
  analyses.index(na) = 1;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if CREATE_MEAN_EPI
  disp('Will create mean resliced EPI image ...')
  na = na+1;
  analyses.type(na) = strmatch('realign', anal_list,'exact');
  analyses.index(na) = 2;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if COREG_MEAN_EPI
  disp('Will coregister mean EPI with coplanars ...')
  na = na+1;
  analyses.type(na) = strmatch('coreg', anal_list,'exact');
  analyses.index(na) = 1;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if COREG_FIRST_EPI
  disp('Will coregister first EPI of first run with coplanars ...')
  na = na+1;
  analyses.type(na) = strmatch('coreg', anal_list,'exact');
  analyses.index(na) = 2;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if COREG_COPLANAR_2_HIRES
   disp('Will coregister the coplanars to the hires  ...')
   na = na+1;
   analyses.type(na) = strmatch('coreg', anal_list,'exact');
   analyses.index(na) = 3;
   analyses.work_dir(na) = 1;
   analyses.mfile(na) = 2;
end  

if COREG_HIRES_2_TEMPLATE
   disp('Will coregister the hires to the template ...')
   na = na+1;
   analyses.type(na) = strmatch('coreg', anal_list,'exact');
   analyses.index(na) = 4; 
   analyses.work_dir(na) = 1;
   analyses.mfile(na) = 2;
end
   
if RESLICE_EPI
  disp('Will reslice EPI images ...')
  na = na+1;
  analyses.type(na) = strmatch('realign', anal_list, 'exact');
  analyses.index(na) = 3;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end


if NORM_COPLANAR_2_TEMPLATE
  disp('will normalize coplanars to template (parameter&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 2;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if NORM_EPI_2_TEMPLATE
  disp('will normalize EPIs to template (parameters&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 3;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end

if NORM_EPI_2_TEMPLATE_CUST_BB
  disp('will normalize EPIs to template using custom bounding box(parameters&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 4;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;  
end

if NORM_EPI_2_TEMPLATE_SPM99_DEF_BB
  disp('will normalize EPIs to template using SPM99 default bounding box (parameters&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 5;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;  
end

if NORM_EPI_2_TEMPLATE_SPM99_BB4
  disp('will normalize EPIs to template using SPM99 bounding box 4 (parameters&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 6;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;  
end

if NORM_HIRES_2_TEMPLATE 
  disp('will normalize hires to template (parameter&reslice) ... ')
  na = na+1;
  analyses.type(na) = strmatch('normalize', anal_list, 'exact');
  analyses.index(na) = 1;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 2;
end
  
if SMOOTH_rEPI
  disp('will smooth resliced EPIs (individuals) ... ')
  na = na+1;
  analyses.type(na) = strmatch('smooth', anal_list, 'exact');
  analyses.index(na) = 1;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 3;
end

if SMOOTH_nEPI 
  disp('will smooth normalized EPIs ... ')
  na = na+1;
  analyses.type(na) = strmatch('smooth', anal_list, 'exact');
  analyses.index(na) = 2;
  analyses.work_dir(na) = 1;
  analyses.mfile(na) = 3;
end

if RUN_MODEL_INDIV
  disp('will run model - individual data ...')

  if SET_F_EQUAL1
    disp('Will set F=1 for storing data to Y.mad ...')
    na = na+1;
    analyses.type(na) = strmatch('defaults_edit', anal_list,'exact');
    analyses.index(na) = 2;
    analyses.work_dir(na) = 1;
    analyses.mfile(na) = 1;   
  end
  
  %use_model_proto = 1;
  nsub_proc = length(proc_subs);
  
  for isub = 1:nsub_proc
    na = na+1;
    analyses.type(na) = strmatch('model', anal_list, 'exact');
    analyses.index(na) = proc_subs(isub);
    work_dir_idx = 2*length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 4;

    % Create analysis directory for single subject data if necessary
    if ~exist(work_dir{work_dir_idx},'dir')
      disp(sprintf('Creating directory: %s', work_dir{work_dir_idx}))
      unix(['mkdir ' work_dir{work_dir_idx}]);
    end
  end
end

if RUN_MODEL_NORM
  disp('will run model - normalized data...')

  if SET_F_EQUAL1
    disp('Will set F=1 for storing data to Y.mad ...')
    na = na+1;
    analyses.type(na) = strmatch('defaults_edit', anal_list,'exact');
    analyses.index(na) = 2;
    analyses.work_dir(na) = 1;
    analyses.mfile(na) = 1;   
  end
  
  %use_model_proto = 1;  
  nsub_proc = length(proc_subs)
  
  for isub = 1:nsub_proc
    na = na+1;
    analyses.type(na) = strmatch('model', anal_list, 'exact');
    analyses.index(na) = proc_subs(isub);
    work_dir_idx = length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 5;

    % Make sure we have a subject directory
    sub_dir_name = fullfile(rootpath, subject_ids{proc_subs(isub)});
    if ~exist(sub_dir_name,'dir')
      disp(sprintf('Creating directory: %s', sub_dir_name))
      unix(['mkdir ' sub_dir_name]);
    end
    
    % Create normalized directory if necessary
    if ~exist(work_dir{work_dir_idx},'dir')
      disp(sprintf('Creating directory: %s', work_dir{work_dir_idx}))
      unix(['mkdir ' work_dir{work_dir_idx}]);
    end
  end
end

if COMPUTE_CONTRASTS_INDIV
  disp('will compute contrasts - individual data ...')

  nsub_proc = length(proc_subs);
  
  for isub = 1:nsub_proc
    na = na+1;
    analyses.type(na) = strmatch('contrasts', anal_list, 'exact');
    analyses.index(na) = proc_subs(isub);
    work_dir_idx = 2*length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 6;                 
%   analyses.mfile(na) = 6+use_model_proto-1;
  end

end

if COMPUTE_CONTRASTS_NORM
  disp('will compute contrasts - normalized data ...')

  nsub_proc = length(proc_subs);
  
  for isub = 1:nsub_proc
    na = na+1;
    analyses.type(na) = strmatch('contrasts', anal_list, 'exact');
    analyses.index(na) = proc_subs(isub);
    work_dir_idx = length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 6;                 
%   analyses.mfile(na) = 6+use_model_proto-1;
  end

end

% NOTE: COMPUTE_GROUP has to come after models and contrasts otherwise working
% directories get screwed up
if COMPUTE_GROUP
    disp('will compute group statistics - normalized data ...')
    ngroup_anal = length(group_analyses);
  
  for ig = 1:ngroup_anal
    
    % Set up the stats
    na = na+1;
    analyses.type(na) = strmatch('stats', anal_list, 'exact');
    analyses.index(na) = ig;
    work_dir{end+1} = fullfile(rootpath, 'group_deriv',group_analyses{ig});
    work_dir_idx = length(work_dir);
    analyses.work_dir(na) = work_dir_idx;
    if ~exist(work_dir{work_dir_idx},'dir')
      disp(sprintf('Creating directory: %s', work_dir{work_dir_idx}))
      unix(['mkdir ' work_dir{work_dir_idx}]);
    end
    analyses.mfile(na) = 10;
    
    % Compute the associated contrasts
    na = na+1;
    analyses.type(na) = strmatch('contrasts', anal_list, 'exact');
    analyses.index(na) = nsub+ig;
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 6;
  end
end


if THRESH_SUM
    disp('will sum thresholded images ...')
    na = na+1;
    analyses.type(na) = strmatch('threshsum', anal_list, 'exact');
    analyses.index(na) = THRESHSUM_CONTRAST;
    analyses.work_dir(na) = 2;
    analyses.mfile(na) = 9;
end

%
% SPECIFY MFILES ASSOCIATED WITH EACH ANALYSIS TYPE
%
% file names are either full pathname or relative to work_dir

mfile = ...
    { fullfile(batch_root, 'simulate_defaults.m'), ...
      fullfile(batch_root, 'simulate_realign.m'), ...      
      fullfile(batch_root, 'simulate_smooth.m'), ...      
      fullfile(batch_root, 'simulate_model_indiv.m'), ...      
      fullfile(batch_root, 'simulate_model_protos.m'), ...      
      fullfile(batch_root, 'simulate_contrasts.m'), ...   
      fullfile(batch_root, 'simulate_group.m'), ...      
      fullfile(batch_root, 'simulate_display.m'), ...      
      fullfile(batch_root, 'simulate_threshsum.m'), ... 
      fullfile(batch_root, 'simulate_stats.m') ... 
      };

