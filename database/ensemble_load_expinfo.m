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
% 2010.06.14 PJ - call to mysql_get_subinfo now requirese that encryption
%                 key information be passed along
% 2010.06.18 FB - now accepts inData={'task','return_outdir'}, for
%                   ensemble_jobman_parallel
% 2011.11.23 PJ - added filtering of subject info

if( isstr( indata ) && strcmp( indata, 'getDefaultParams' ) )
	result.ensemble.experiment_title = 'use local settings';
	result.mysql.host = 'use local settings';
	result.mysql.database = 'use local settings';
	result.mysql.conn_id = -1;
	return;
end

if (iscell(indata) && ~isempty(indata) && isfield(indata{1},'task') && ...
        ~isempty(strmatch('return_outdir',indata{1}.task))) || ...
        (isstruct(indata) && isfield(indata,'task') && ...
        ~isempty(strmatch('return_outdir',indata.task)))
    result = '';
    return
end

% Initialize the output data struct
result = ensemble_init_data_struct;
result.type = 'experiment_info';

% Check to see if we have a valid connection
if ~isfield(params,'mysql')
	params.mysql = params.ensemble;
end
conn_id = mysql_make_conn(params.mysql);

expname = params.ensemble.experiment_title;

%% Get the experiment metadata
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

%% Load session table info
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

%% Load response table
fprintf('Loading response table information ...\n');
respinfo = ensemble_init_data_struct;
[respinfo.data,respinfo.vars] = mysql_extract_data('table', exp_meta.response_table, ...
	'extract_flds', extract_vars, 'order_by','response_id', ...
	'session_id', sessIDList, ...
	'conn_id', conn_id);
respinfo.name = 'response_data';
respinfo.type = 'response_data';
respcols = set_var_col_const(respinfo.vars);

% Check to make sure that subquestion information isn't NaN
if any(isnan(respinfo.data{respcols.subquestion}))
	fprintf('Replacing NaNs in subquestion field with 1\n');
	respinfo.data{respcols.subquestion}(isnan(respinfo.data{respcols.subquestion})) = 1;
end

% Check to see if we are filtering using 'terminal_form'
respFilt = [];
try term_form = params.ensemble.terminal_form; catch term_form = []; end
if ~isempty(term_form) && remove_incomplete
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
if remove_incomplete && ~isempty(term_form)
	fprintf('Removing incomplete sessions from session_info ...\n');
	clear tmpFilt
	tmpFilt.include.any.session_id = finished_sessids;
	sessInfo = ensemble_filter(sessInfo, tmpFilt);
end

% update the output structure
result.vars{end+1} = 'session_info';
result.data{end+1} = sessInfo;
outCols = set_var_col_const(result.vars);

%% Load subject table info
% Pull the subject information for all of the subjects
if isfield(params, 'ignoreSubjectTable') && params.ignoreSubjectTable
	fprintf('%s: ignoring subject table information\n', mfilename);
else
	fprintf('Loading subject information ...\n');
	if ~isfield(params.ensemble,'enc_key')
		if isfield(params, 'mysql') && isfield(params.mysql, 'enc_key')
			params.ensemble.enc_key = params.mysql.enc_key;
		else
			params.ensemble.enc_key = '';
		end
	end
	
	if ~isfield(params.ensemble, 'conn_id') || isempty(params.ensemble.conn_id)
		try
			params.ensemble.conn_id = params.mysql.conn_id;
		catch
			fprintf('Failed to locate connection ID');
		end
	end
	
	if isfield(params,'mysql')
		login_params = params.mysql;
	elseif isfield(params,'ensemble')
		login_params = params.ensemble;
	else
		login_params = struct();
	end
	
	subInfo = mysql_get_subinfo('subject_id', sessInfo.data{sess_colNames.subject_id}, ...
		'mysql', login_params);
	subInfo.name = 'subject_info';
	subInfo.type = 'subject_info';
	
	% Perform filtering on subject info
	if isfield(params,'filt')
		subInfo = ensemble_filter(subInfo, params.filt);
	end
	scols = set_var_col_const(subInfo.vars);
	
	% Get list of retained subject IDs
	goodSubIDs = unique(subInfo.data{scols.subject_id});
	
	% Refilter session and response table structures
	if isfield(params,'ensemble') && isfield(params.ensemble,'refilterBySubjectTableIDs') && ~params.ensemble.refilterBySubjectTableIDs
		fprintf('Not refiltering session and response data using subject table information ...\n');
	else
		fprintf('Refiltering session and response data using subject table information ...\n');
		tmpfilt = struct();
		tmpfilt.include.all.subject_id = goodSubIDs;
		
		result.data{outCols.session_info} = ...
			ensemble_filter(result.data{outCols.session_info}, tmpfilt);
		sessInfo = result.data{outCols.session_info};
		
		result.data{outCols.response_data} = ...
			ensemble_filter(result.data{outCols.response_data}, tmpfilt);
		respinfo = result.data{outCols.response_data};
	end
	
	% update the output structure
	result.vars{end+1} = 'subject_info';
	result.data{end+1} = subInfo;
end % if ignoreSubjectTable

%% Get stimulus info
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


