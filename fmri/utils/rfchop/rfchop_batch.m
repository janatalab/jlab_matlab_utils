% 
%  Specifications for batch processing data from selective attention to
%  timbre study (ts1).
%

% 5/15/00 PJ - coded for Schubert study
% 6/16/00 PJ - modfied for ts1 study
% 2/19/01 PJ - modified for re-analysis.  Among other things, the model now has
% a peak at the beginning of the attend conditions to better capture the
% attentional orienting phase that most subjects report experiencing at the
% beginning of the excerpt

%
%  GLOBAL variables
%

global subject_ids proc_subs use_model_proto sinfo
global rootpath

global COREG_WITH_MUTUAL_INFO		% Perform coregistration using mutual
                                        % information method rather than the
                                        % standard method.  This hopefully
                                        % overcomes some of the distortions due
                                        % to the FIL EPI template.

global model_action regressors2add created_model_specs % these are initialized by RUN_MODEL
global group_analyses			% group statistics analyses

rootpath = '/data1/ts1/';

%
%  STRUCTURE CONTAINING LIST OF ANALYSES TO PERFORM
%
%  The analyses structure is created dynamically based on the list of desired
%  analyses
%

SET_GENERAL_DEFS = 1;

%
% PREPROCESSING OPTIONS
%

REALIGN_EPI = 0;
CREATE_MEAN_EPI = 0;

COREG_WITH_MUTUAL_INFO = 1;
COREG_MEAN_EPI = 0;
COREG_FIRST_EPI = 0;			% I don't usually do this

COREG_COPLANAR_2_HIRES = 0;
COREG_HIRES_2_TEMPLATE = 0;			% I don't usually do this

RESLICE_EPI = 0;			% I don't usually do this

NORM_COPLANAR_2_TEMPLATE = 0;			% I don't usually do this

NORM_EPI_2_TEMPLATE = 0;
NORM_EPI_2_TEMPLATE_CUST_BB = 0;
NORM_EPI_2_TEMPLATE_SPM99_DEF_BB = 0;	% This is the one I usually choose
NORM_EPI_2_TEMPLATE_SPM99_BB4 = 0;

NORM_HIRES_2_TEMPLATE = 0;	% You have to do this if your voxel size in
                                % normalized hires volumes is smaller than the
                                % voxel size in normalized EPI volumes

SMOOTH_rEPI = 0;   % I don't usually do this because I haven't resliced
SMOOTH_nEPI = 0;

%
% MODELLING OPTIONS
%

RUN_MODEL_INDIV = 0;			% Create and evaluate model of single
                                        % subject, unnormalized data
RUN_MODEL_NORM = 0;			% Create and evaluate model of single
                                        % subject, normalized data (for group analyses)
	
model_action = 'spec_and_estimate'; %{'specify','review','estimate','spec_and_estimate'}
regressors2add = {'linear','motion'};
					
COMPUTE_CONTRASTS_INDIV = 0;
COMPUTE_CONTRASTS_NORM = 1;

COMPUTE_GROUP = 1;
group_outdir = 'group';
group_analysis_list = {'_all_t'};	% specify _all_t to do all

THRESH_SUM = 0; THRESHSUM_CONTRAST=2;

%
%  SUBJECT INFO
%
sinfo = get_ts1_sinfo;

subject_ids = cellstr(char(sinfo(:).id))

nsub = length(subject_ids);
proc_subs = [1:12]			% 1:12 - original, 13:17 - thin slices

%proc_subs = [14 15];				% try non mutual info coreg w/IO and RP

%
%  LIST OF ANALYSIS TYPES
%

type = {'defaults_edit','model','contrasts','headers',...
'means','realign','coreg','normalize','smooth','display','threshsum','stats'};

%
%  GROUP ANALYSIS PARAMETERS
%

anal_list = strvcat(type);

analyses = struct('type',[],'index',[],'work_dir',[],'mfile',[]);

na = length(analyses.type);

%
%  LIST OF WORK DIRECTORIES
%

work_dir = {rootpath, fullfile(rootpath, 'group'), ...
      };

for isub = 1:nsub
  work_dir{end+1} = fullfile(rootpath, subject_ids{isub}, 'single');
end
 
for isub = 1:nsub
  work_dir{end+1} = fullfile(rootpath, subject_ids{isub}, 'normalized');
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

  use_model_proto = 1;
  
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

  use_model_proto = 2;

  created_model_specs = 0;
  
  nsub_proc = length(proc_subs);
  
  for isub = 1:nsub_proc
    na = na+1;
    analyses.type(na) = strmatch('model', anal_list, 'exact');
    analyses.index(na) = proc_subs(isub);
    work_dir_idx = length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 5;

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
    analyses.index(na) = 1;
    work_dir_idx = 2*length(work_dir) - nsub + proc_subs(isub);
    analyses.work_dir(na) = work_dir_idx;
    analyses.mfile(na) = 6;
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
  end

end

% NOTE: COMPUTE_GROUP has to come after models and contrasts otherwise working
% directories get screwed up
if COMPUTE_GROUP
  disp('will compute group statistics - normalized data ...')
    
  % Make sure the group directory exists
  groupdir = fullfile(rootpath, group_outdir);
  if ~exist(groupdir,'dir')
    disp(sprintf('Creating directory: %s', groupdir))
    unix(['mkdir ' groupdir]);
  end

  % if group_analyses is specified as '_all_t', then execute the contrasts
  % script to get all the names of the T contrasts and expand group_analyses to
  % include all of these
  
  ts1_contrasts
  
  if strmatch('_all_t', group_analysis_list, 'exact')
    ng = 0;
    for ic = 1:length(contrasts(1).names)
      if contrasts(1).types{ic} == 'T'
	ng = ng+1;
	group_analyses{ng} = contrasts(1).names{ic};
      end
    end
  else
    group_analyses = group_analysis_list;
  end
  
  ngroup_anal = length(group_analyses);
  
  for ig = 1:ngroup_anal
    
    % Set up the stats
    na = na+1;
    analyses.type(na) = strmatch('stats', anal_list, 'exact');
    analyses.index(na) = ig;
    work_dir{end+1} = fullfile(groupdir,group_analyses{ig});
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
    { '/data1/matlab/ts1/ts1_defaults.m', ...
      '/data1/matlab/ts1/ts1_realign.m', ...      
      '/data1/matlab/ts1/ts1_smooth.m', ...      
      '/data1/matlab/ts1/ts1_model_indiv.m', ...      
      '/data1/matlab/ts1/ts1_model_norm.m', ...      
      '/data1/matlab/ts1/ts1_contrasts.m', ...      
      '/data1/matlab/ts1/ts1_group.m', ...      
      '/data1/matlab/ts1/ts1_display.m', ...      
      '/data1/matlab/ts1/ts1_threshsum.m', ...      
      '/data1/matlab/ts1/ts1_stats.m' ...
      };

