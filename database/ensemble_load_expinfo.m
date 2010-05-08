function result = ensemble_load_expinfo(indata,params) 
% returns the response table and metadata associated with an experiment
%
%   result = ensemble_load_expinfo(indata,params)
% 
% REQUIRES
%   params.ensemble.conn_id
%   params.ensemble.experiment_title
%   params.ensemble.remove_incomplete_sessions
%   params.ensemble.terminal_form
%   params.extract_vars
%   params.filt
% 
% the 'indata' variable is not used yet, but is simply there as a
% placeholder. It may be used in future versions of this function
%

% Original version by Petr Janata
% 16 Feb 2007 - reorganized code by PJ into a function, Stefan Tomic
%
% 15 March 2007 - edited function to conform to ensemble data structures
% 04 April 2007 - PJ implemented incremental buildup of result structure.
%                 Added check to make sure there are stimuli to load
%                 information about.
% 03 May 2007 - Added option to pass in a list of variables to extract from the
%               response table.
% 05/18/07 PJ - Added handling of terminal_form which is used in conjunction
%               with remove_incomplete_sessions to specify experiment
%               completion relative to a particular form ID rather than an
%               end_datetime entry in the session table.
%
% 08/06/08 JG - Added getDefaultParams feature
%
% 2009.05.18 FB - added some additional header documentation
% 2009.10.18 PJ - added check to make sure response table was found
% 2010.01.20 FB - now also returns 'misc_info' in the 'response_data' struct
% 2010.05.08 PJ - isolated subject_summary_stats in try/catch

if( isstr( indata ) && strcmp( indata, 'getDefaultParams' ) )
	result.ensemble.experiment_title = 'use local settings';
	result.ensemble.host = 'use local settings';
	result.ensemble.database = 'use local settings';
	result.ensemble.conn_id = -1;
	return;
end


% Initialize the output data struct
result = ensemble_init_data_struct;

try
  conn_id = params.ensemble.conn_id;
catch
  conn_id = 1;
  tmp_conn_id = 1;
end

expname = params.ensemble.experiment_title;

% Get the experiment information
fprintf('Getting experiment information for: %s\n', expname);
exp_meta = mysql_extract_metadata('table','experiment', ...
  'experiment_title',expname, ...
  'conn_id', conn_id);

if isempty(exp_meta.response_table)
  warning(['Failed to locate response_table for experiment (%s)\n' ...
    'Check to make sure that you are connecting to the database ' ...
    'you think you are connecting to!\n'], expname);
  return
end
        
% Pull desired parts of the response table into a response structure
fprintf('Extracting the response table: %s\n', exp_meta.response_table);
if ~isfield(params.ensemble, 'extract_vars') || isempty(params.ensemble.extract_vars)
  extract_vars = {'session_id','subject_id', ...
	'response_order','response_id', 'date_time', ...
	'form_id','form_order','question_id','subquestion', ...
	'stimulus_id','trial_id','response_enum','response_text','misc_info'};
else
  extract_vars = params.ensemble.extract_vars;
end

exp_meta_dataStruct = ensemble_tree2datastruct(exp_meta);
exp_meta_dataStruct.name = 'experiment_metadata';
exp_meta_dataStruct.type = 'experiment_metadata';

% update the output structure
result.vars{end+1} = 'experiment_metadata';
result.data{end+1} = exp_meta_dataStruct;

%get the experiment ID
exp_meta_cols = set_var_col_const(exp_meta_dataStruct.vars);
experiment_id = exp_meta_dataStruct.data{exp_meta_cols.experiment_id};

%get sessinfo data
sessInfo = ensemble_init_data_struct;
sessInfo.name = 'session_info';
sessInfo.type = 'session_info';
[sessInfo.data,sessInfo.vars] = mysql_extract_data('table','session', ...
    'experiment_id',experiment_id, ...
    'conn_id',conn_id);

%apply any filters to the sessions. This could trap subjects to exclude
if(isfield(params,'filt'))
  sessFilt = params.filt;
  sessInfo = ensemble_filter(sessInfo,sessFilt);
end

%filter out any session that have NaN as the end_datetime (these
%sessions were not completed)
remove_incomplete = 0;
if(isfield(params.ensemble,'remove_incomplete_sessions') && ...
      params.ensemble.remove_incomplete_sessions == 1)
  remove_incomplete = 1;
  if ~isfield(params.ensemble,'terminal_form') || isempty(params.ensemble.terminal_form)
    sessFilt.exclude.any.end_datetime = NaN;
    sessInfo = ensemble_filter(sessInfo,sessFilt);
  end
end

%get the filtered list of session IDs
sess_colNames = set_var_col_const(sessInfo.vars);
sessIDList = sessInfo.data{sess_colNames.session_id};

%get the response info and filter it
fprintf('Loading response table information ...\n');
respinfo = ensemble_init_data_struct;
[respinfo.data,respinfo.vars] = mysql_extract_data('table', exp_meta.response_table, ...
    'extract_flds', extract_vars, 'order_by','response_id', ...
    'session_id', sessIDList, ...
    'conn_id', conn_id);
respinfo.name = 'response_data';
respinfo.type = 'response_data';
respcols = set_var_col_const(respinfo.vars);

% Check to see if we are filtering using 'terminal_form'
respFilt = [];
try term_form = params.ensemble.terminal_form; catch term_form = []; end
if ~isempty(term_form) & remove_incomplete
  form_mask = respinfo.data{respcols.form_id} == term_form;
  finished_sessids = unique(respinfo.data{respcols.session_id}(form_mask));
  respFilt.include.any.session_id = finished_sessids;
end

if isfield(params,'filt')
  respFilt = add2filt(respFilt,params.filt);
end

if(~isempty(respFilt))
  fprintf('Filtering response table ...\n');
  respinfo = ensemble_filter(respinfo,respFilt);
end

% update the output structure
result.vars{end+1} = 'response_data';
result.data{end+1} = respinfo;

% Refilter session information if necessary
if remove_incomplete & ~isempty(term_form)
  fprintf('Removing incomplete sessions from session_info ...\n');
  clear tmpFilt
  tmpFilt.include.any.session_id = finished_sessids;
  sessInfo = ensemble_filter(sessInfo, tmpFilt);
end

% update the output structure
result.vars{end+1} = 'session_info';
result.data{end+1} = sessInfo;

% Pull the subject information for all of the subjects
fprintf('Loading subject information ...\n');
subInfo = mysql_get_subinfo('subject_id', ...
			    sessInfo.data{sess_colNames.subject_id},'conn_id', conn_id);
subInfo.name = 'subject_info';
subInfo.type = 'subject_info';

% update the output structure
result.vars{end+1} = 'subject_info';
result.data{end+1} = subInfo;

% Get stimulus information if there are any stimuli associated with the
% response table
if ~all(isnan(respinfo.data{respcols.stimulus_id}))
    fprintf('Loading stimulus and attribute information...\n');
    stimulusDataStruct = ensemble_get_stiminfo(respinfo, params);
    result.vars = [result.vars stimulusDataStruct.vars];
    result.data = [result.data stimulusDataStruct.data];
end 

%get subject summary stats
try
  subSummaryStats = ensemble_summary_subject_stats({subInfo sessInfo});
  result.vars{end+1} = 'subject_summary_stats';
  result.data{end+1} = subSummaryStats;
catch
  warning('Unable to run ensemble_summary_subject_stats');
end

if(exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
