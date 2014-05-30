function outData = check_completion_v2(varargin)
% provides details on subjects who have started and who have completed a given experiment
%
% adapted from anps_collect_check_completion.m written by Fred Barrett
%
% REQUIRED PARAMS
% Either 
%    params.resptbl_name              - The name of the response table for the experiment
% or 
%    params.ensemble.experiment_title - The title of the experiment
%
% OPTIONAL PARAMS
% params.terminal_form             - The form ID to check against for
%                                    completion of experiment
% params.terminal_question         - The question ID to check against for
%                                    completion
% params.filt                      - Filter params to filter out unwanted subjects
% params.ensemble.conn_id          - Connection ID to use for database 
%                                    if no live connection, a new connection is established
% params.ensemble.host             - hostname of database
% params.ensemble.database         - database name
% params.report_after              - any sessions before this date will be
%                                    ignored in the report (if not specified,
%                                    all sessions will display)
% params.check_completion_v2.report.incomplete_details - more detailed
%                                    information about incomplete sessions
%
% OUTPUTS
%
% Returns a list of subjects who completed the experiment

% Author(s)
% 2/21/2009 Fred Barrett - Original Author
% 4/7/2009 Stefan Tomic - adapted script into a more general purpose function for all experiments
% 06/15/2010 PJ - eliminated ensemble_globals, and forced new
%                 msyql_make_conn scheme.
% 13Mar2012 PJ - Fixed handling of database connectivity, giving precedence to
%                mysql. Added reporting of subject names and start and stop
%                times
% 03Nov2013 PJ - Added optional reporting of more details for incomplete sessions

outData = [];

if nargin == 1
  params = varargin{1};
elseif nargin == 2
  params = varargin{2};
else
  fprintf('%s: Wrong number of arguments: %d\n', mfilename, nargin);
end

if isfield(params,'mysql')
	db_params_struct = 'mysql';
elseif isfield(params, 'ensemble')
	db_params_struct = 'ensemble';
else
	error('Do not have a mysql structure with db parameter info');
end

try
  conn_id = params.mysql.conn_id;
catch
	try 
		conn_id = params.ensemble.conn_id;
	catch
		conn_id = 7;
		params.mysql.conn_id = conn_id;
	end
end

try
  params.filt;
catch
  params.filt = {};handling
end

if(mysql(conn_id,'status') ~= 0)
  mysql_make_conn(params.(db_params_struct));
end

try tbl_name = params.resptbl_name; catch tbl_name = ''; end

if isempty(tbl_name) && isfield(params.ensemble, 'experiment_title')
  mysql_str = sprintf('SELECT response_table, experiment_id FROM experiment WHERE experiment_title="%s"', params.ensemble.experiment_title);
  [tbl_name, exp_id] = mysql(conn_id, mysql_str);
  tbl_name = tbl_name{1};
  if ~isfield(params.ensemble, 'response_table')
    params.ensemble.response_table = tbl_name;
  end
end

try form_id = params.terminal_form; catch form_id = []; end;

% If form_id is unspecified, get it from the experiment_x_form table
if isempty(form_id)
  mysql_str = sprintf('SELECT form_id, form_order FROM experiment_x_form WHERE experiment_id=%d;', exp_id);
  [ids,order] = mysql(conn_id, mysql_str);
  % Get next to last form, as last form is an end_session form that doesn't
  % leave its mark on the response table
  form_id = ids(order==(max(order)-1));
end

try question_id = params.terminal_question; catch question_id = []; end
if isempty(question_id)
  mysql_str = sprintf('SELECT question_id, form_question_num FROM form_x_question WHERE form_id=%d;', form_id);
  [qids, order] = mysql(conn_id,mysql_str);
  question_id = qids(order==max(order));
end

try PRINT_TO_FILE = params.report.write2file; catch PRINT_TO_FILE=0; end

% % % REPORT AFTER THIS DATE
try
  report_after = datenum(params.report_after);
catch
  %if report_after wasn't specied, report all sessions
  report_after = 0;
end

%% Consult the session table to get list of completed subs, sessions, and time stamps
vars = {'session_id','date_time','end_datetime','subject_id','ticket_id'};
mysql_str = sprintf('SELECT %s FROM session WHERE experiment_id =%d;', cell2str(vars,','), exp_id);
data = cell(1,length(vars));
[data{:}] = mysql(conn_id,mysql_str);
completedData = ensemble_init_data_struct;
completedData.type = 'session_info';
completedData.vars = vars;
completedData.data = data;
cdCols = set_var_col_const(vars);

% Get the ticket codes

%% Filter completion data if desired
completedData = ensemble_filter(completedData,params.filt);

% open logfile
if exist('PRINT_TO_FILE','var') && PRINT_TO_FILE
  logfname = fullfile(params.paths.logpath,'completion_info.txt');
  fid = fopen(logfname,'wt');
  fprintf(fid,'Completion information for experiment: %s\n', params.ensemble.experiment_title);
  fprintf(fid,'Generated: %s\n\n\n', datestr(now));
else
  fid = 1; % stdout
end
params.report.fid = fid;

% Get subject information
subInfo = mysql_get_subinfo('subject_id', completedData.data{cdCols.subject_id}, ...
  'mysql', params.(db_params_struct));
siCols = set_var_col_const(subInfo.vars);

nsess = length(completedData.data{cdCols.session_id});
fprintf(fid,'%d initiated sessions\n', nsess);

completed_mask = ~isnan(completedData.data{cdCols.end_datetime});
completed_idxs = find(completed_mask);
fprintf(fid,'%d completed sessions\n', sum(completed_mask));

% Add the completion mask to the completedData structure
completedData.vars{end+1} = 'completed';
completedData.data{end+1} = completed_mask;
cdCols = set_var_col_const(completedData.vars);

for itype = 1:2
  if itype == 1
    idx_list = completed_idxs;
    type_str = 'COMPLETED';
  else
    idx_list = find(~completed_mask);
    type_str = 'INCOMPLETE';
  end
  
  fprintf(fid,'\n%s: %d sessions\n', type_str, length(idx_list));
  fprintf(fid,'Session\tSubject\tFirstName\tLastName\tStartTime\tEndTime\tTicketID\tTicketCode\n');
  for isess = 1:length(idx_list)
    curr_idx = idx_list(isess);
    
    curr_tick_id = completedData.data{cdCols.ticket_id}(curr_idx);
    mysql_str = sprintf('SELECT ticket_code FROM ticket WHERE ticket_id = %d;', curr_tick_id);
    ticket_code = mysql(conn_id, mysql_str);
    
		start_str = datestr(completedData.data{cdCols.date_time}(curr_idx));
		if strcmp(type_str,'INCOMPLETE')
			stop_str = '';
		else
			stop_str = datestr(completedData.data{cdCols.end_datetime}(curr_idx));
		end
		
		first_name = subInfo.data{siCols.name_first}{strcmp(subInfo.data{siCols.subject_id},completedData.data{cdCols.subject_id}{curr_idx})};
		last_name = subInfo.data{siCols.name_last}{strcmp(subInfo.data{siCols.subject_id},completedData.data{cdCols.subject_id}{curr_idx})};
			
    fprintf(fid, '%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n', ...
      completedData.data{cdCols.session_id}(curr_idx), ...
      completedData.data{cdCols.subject_id}{curr_idx}, ...
			first_name, last_name, ...
      start_str, stop_str, ...
      curr_tick_id, ticket_code{1});
  end % for isess
end % for itype (COMPLETED, INCOMPLETE)

%
% Check to see if there is any further reporting to be done
%
if isfield(params, mfilename) && isfield(params.(mfilename), 'report')
  report_types = fieldnames(params.(mfilename).report);
  nreports = length(report_types);
  for ireport = 1:nreports
    currReport = report_types{ireport};
    fprintf('\n\nRunning report %d/%d: %s\n', ireport, nreports, currReport);
    
    % Get a handle to the sub-function we want to run
    fh = str2func(currReport);
    result = fh({completedData, subInfo}, params);
  end
end % if isfield(params, mfilename) && isfield(params.(mfilename), 'report')

return
end

function result = incomplete_details(data_st, params)
  st = dbstack;
  funcname = st(1).name; % get current function name
  
  % Get the parameters that specifically control this function
  fparams = params.(mfilename).report.(funcname);
  if isnumeric(fparams) && ~fparams
    result = [];
    return
  end
  
  if isstruct(fparams)
    if isfield(fparams, 'fname')
      fid = ensemble_init_fid(fparams);
    end
  else
    fid = 1;
  end
  
  % Parse out the incoming data structures
  an_idx = ensemble_find_analysis_struct(data_st,struct('type','session_info'));
  sess_st = data_st{an_idx};
  sessCols = set_var_col_const(sess_st.vars);
  
  an_idx = ensemble_find_analysis_struct(data_st,struct('type','subject_info'));
  sub_st = data_st{an_idx};
  subCols = set_var_col_const(sub_st.vars);
  
  % Filter out the completed sessions
  cfilt.exclude.all.completed = 1;
  sess_st = ensemble_filter(sess_st, cfilt);
  sessIDs = sess_st.data{sessCols.session_id};
  nsess = length(sessIDs);
  
  % Get the response table data for the incomplete sessions
  respdata = ensemble_init_data_struct;
  [respdata.data, respdata.vars] = mysql_extract_data(...
    'table', params.ensemble.response_table, ...
    'conn_id', params.mysql.conn_id, ...
    'session_id', sess_st.data{sessCols.session_id});
  respCols = set_var_col_const(respdata.vars);
  
  % Loop over sessions
  lastStim = nan(nsess,1);
  for isess = 1:nsess
    currSessID = sessIDs(isess); 
    currSubID = sess_st.data{sessCols.subject_id}{isess};
    
    % Get mask into subject info structure
    submask = strcmp(sub_st.data{subCols.subject_id},currSubID);
    
    % Filter the response table data to extract the current session
    cfilt = [];
    cfilt.include.all.session_id = currSessID;
    currData = ensemble_filter(respdata,cfilt);
    
    % Get the last entry that was made for this subject
    lastTimeStamp = max(currData.data{respCols.date_time});
    lastEntryMask = currData.data{respCols.date_time} == lastTimeStamp;
    
    % Print out diagnostic information
    fprintf(fid,'\nSession ID: %d, Subject ID: %s, Name: %s, Started: %s\n', ...
      currSessID, ...
      currSubID, ...
      sprintf('%s %s', sub_st.data{subCols.name_first}{submask}, sub_st.data{subCols.name_last}{submask}), ...
      datestr(sess_st.data{sessCols.date_time}(isess)));
    
    % Provide information about last form submitted
    formID = currData.data{respCols.form_id}(lastEntryMask);
    if any(diff(formID))
      error('Multiple formIDs associated with single timestamp: %s', sprintf('\t%d', unique(formID)))
    else
      formID = formID(1);
    end
    
    mysql_str = sprintf('SELECT form_name FROM form WHERE form_id = %d;', formID);
    formName = mysql(params.mysql.conn_id, mysql_str);
    formName = formName{1};
    
    fprintf(fid,'\tLast submitted form: ID=%d, %s, %s\n', formID, formName, datestr(lastTimeStamp));
    
    % Provide number and list of stimuli encountered
    stimIDs = unique(currData.data{respCols.stimulus_id}, 'stable');
    stimIDs(isnan(stimIDs)) = [];  % eliminate entries for which there were no stimulus IDs
    fprintf(fid,'\t%d stimuli:%s\n', length(stimIDs), sprintf('\t%d', stimIDs));
    
    % Track last stimulus encountered across sessions 
    if ~isempty(stimIDs)
      lastStim(isess) = stimIDs(end);
    end
    
  end % for isess
  result = 0;
end % incomplete_details

