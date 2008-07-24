function [data] = extract_generic(resp_tbl,form_id,extract_vars, subids, conn_id)
% [data] = extract_generic(resp_tbl, form_id, extract_vars, subids, conn_id);
%
% Generic data extraction for data associated with a specfic form in an
% experiment response table in the database.
%
% extract_vars - cell string containing optional list of fields to extract
%
% subids - subjects to process
%

% 03/07/05 PJ

global SQL_HOST DATABASE

% Check for connection to database
try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
end

try 
  extract_vars{1}; 
catch
  extract_vars = {...
	'subject_id', ...
	'stimulus_id', ...
	'date_time', ...
	'question_id', ...
	'subquestion', ...
	'question_iteration', ...
	'response_id', ...
	'response_enum', ...
	'response_text' ...
    };
end
    
data.vars = extract_vars;

try subids(1);
  subject_str = sprintf('subject_id="%s" OR ', subids{:});
  subject_str(end-3:end) = [];
  subject_str = sprintf('AND (%s)', subject_str);
catch
  subject_str = '';
end
  
% Extract the data
sql_str = sprintf('select %s from %s where form_id=%d %s ;', cell2str(extract_vars,','), resp_tbl, form_id, subject_str);
[data.data{1:length(extract_vars)}] = mysql(conn_id,sql_str);

% Resolve the question IDs to get the associated question text
sql_str = sprintf(['SELECT question_text FROM question,%s WHERE' ...
      ' question.question_id=%s.question_id AND %s.form_id=%d %s ;'], resp_tbl, resp_tbl, resp_tbl, form_id, subject_str);
[qtxt] = mysql(conn_id,sql_str);
data.vars{end+1} = 'question_text';
data.data{end+1} = qtxt;

% Get the subquestion text
sql_str = sprintf(['SELECT heading FROM' ...
      ' question_x_data_format,%s WHERE' ...
      ' question_x_data_format.question_id=%s.question_id AND' ...
      ' question_x_data_format.subquestion=%s.subquestion AND %s.form_id=%d %s;'], ...
    resp_tbl, resp_tbl, resp_tbl, resp_tbl, form_id, subject_str);
[heading] = mysql(conn_id,sql_str);
data.vars{end+1} = 'heading';
data.data{end+1} = heading;

if ~conn_id
  mysql(conn_id,'close');
end

return