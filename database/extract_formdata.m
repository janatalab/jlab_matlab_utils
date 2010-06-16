function [data] = extract_formdata(resp_tbl,form_ids,varargin)
% Extracts data for given forms in the given response table.
% [data] = extract_formdata(resp_tbl, form_ids, varargin);
%
% Generic data extraction for data associated with a specfic list of forms in an
% experiment response table in the database. Data is a cell array with each
% form in one element
%
% A flexible number of options can be passed into the function as tag/value
% pairs.
%
% 'conn_id', conn_id - the ID of the mysql connection - REQUIRED
% 'subids', {subids} - a cell array of subject IDs
% 'sessids', [sessids] - a vector of session IDs
% 'extract_vars', {extract_vars} - cell array of field names to
%                                  extract. Default is all.

% 06/22/06 PJ - Adapted from extract_generic
% 10/28/06 PJ - enabled handling of vector of form_ids. 
% 06/15/10 PJ - sanitized mysql_make_conn

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case 'conn_id'
      conn_id = varargin{iarg+1};
    case 'subids'
      subids = varargin{iarg+1};
    case 'sessids'
      sessids = varargin{iarg+1};
    case 'extract_vars'
      extract_vars = varargin{iarg+1};
    otherwise
      fprintf('extract_formdata: Unknown input argument: %s\n', varargin{iarg});
  end
end

if ~exist('conn_id','var') || isempty(conn_id)
  error('%s: Do not have a valid connection ID', mfilename);
end

try 
  extract_vars{1}; 
catch
  % This should be made more dynamic by getting a list of all field names in
  % the desired response table.
  extract_vars = {...
	'subject_id', ...
	'session_id', ...
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
    
try subids(1);
  subject_str = sprintf('"%s",', subids{:});
  subject_str(end) = [];
  subject_str = sprintf('AND subject_id IN (%s)', subject_str);
catch
  subject_str = 'AND 1';
end
  
try sessids(1);
  sess_str = sprintf('%d,', sessids);
  sess_str(end) = [];
  sess_str = sprintf('AND session_id IN (%s)', sess_str);
catch
  sess_str = 'AND 1';
end

nform = length(form_ids);
data = cell(nform,1);
for iform = 1:nform
  data{iform}.vars = extract_vars;
  form_id = form_ids(iform);
  % Extract the data
  sql_str = sprintf('SELECT %s FROM %s WHERE form_id=%d %s %s;', cell2str(extract_vars,','), resp_tbl, form_id, subject_str, sess_str);
  [tmp{1:length(extract_vars)}] = mysql(conn_id,sql_str);
  data{iform}.data = tmp;

  % Resolve the question IDs to get the associated question text
  sql_str = sprintf(['SELECT question_text FROM question,%s WHERE' ...
	' question.question_id=%s.question_id AND %s.form_id=%d ' ...
	' %s %s;'], resp_tbl, resp_tbl, resp_tbl, form_id, subject_str, sess_str);
  [qtxt] = mysql(conn_id,sql_str);
  data{iform}.vars{end+1} = 'question_text';
  data{iform}.data{end+1} = qtxt;

  % Get the subquestion text
  sql_str = sprintf(['SELECT heading FROM' ...
	' question_x_data_format,%s WHERE' ...
	' question_x_data_format.question_id=%s.question_id AND' ...
	' question_x_data_format.subquestion=%s.subquestion AND %s.form_id=%d %s %s;'], ...
      resp_tbl, resp_tbl, resp_tbl, resp_tbl, form_id, subject_str, sess_str);
  [heading] = mysql(conn_id,sql_str);
  data{iform}.vars{end+1} = 'heading';
  data{iform}.data{end+1} = heading;
end % for iform

if ~conn_id
  mysql(conn_id,'close');
end

return
