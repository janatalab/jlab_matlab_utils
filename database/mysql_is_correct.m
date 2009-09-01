function outstr = mysql_is_correct(varargin)
% Checks to see whether the response to a trial is correct
%
% Input arguments are passed in as parameter/value pairs.  
% They must include:
%
% 'trial_id'
% 'session_id'

% 08/31/2009 Petr Janata

global mysql_conn_id

% Parse input arguments
trial_id = [];
session_id = [];

narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'trial','trial_id'}
      trial_id = varargin{iarg+1};
    case {'session','session_id'}
      session_id = varargin{iarg+1};
    otherwise
      fprintf('%s: Unknown input argument: %s\n', mfilename, varargin{iarg});
  end
end

if isempty(trial_id)
  error('%s: Trial ID not specified', mfilename);
end
if isempty(session_id)
  error('%s: Session ID not specified', mfilename);
end

% Check for connection to database
try mysql_conn_id(1);
  conn_id = mysql_conn_id;
  tmp_conn_id = 0;
catch   
  mysql_make_conn;
  conn_id = 0;
  tmp_conn_id = 1;
end

% We need to get the response table name, based on the session ID
mysql_str = sprintf(['SELECT response_table FROM experiment WHERE ' ...
      'experiment_id = (SELECT experiment_id FROM session WHERE session_id = %d);'], session_id);
response_table = mysql(conn_id, mysql_str);
response_table = response_table{1};

% Get the trial information from the trial table
mysql_str = sprintf(['SELECT correct_response_enum, correct_response_text FROM trial ', ...
      'WHERE trial_id = %d;'], trial_id);
[correct_response_enum, correct_response_text] = mysql(conn_id, mysql_str);
correct_response_text = correct_response_text{1};
if isempty(correct_response_text) && isnumeric(correct_response_enum);
  use_enum = true;
else
  use_enum = false;
end

% Get the response to this session and trial from the response table
mysql_str = sprintf(['SELECT response_enum, response_text FROM %s ' ...
      'WHERE session_id = %d AND trial_id = %d;'], response_table, session_id, trial_id);
[response_enum, response_text] = mysql(conn_id, mysql_str);
response_text = response_text{1};

if use_enum
  if response_enum == correct_response_enum
    outstr = 'TRUE';
  else
    outstr = 'FALSE';
  end
else
  if strcmp(response_text, correct_response_text)
    outstr = 'TRUE';
  else
    outstr = 'FALSE';
  end
end
  
if tmp_conn_id
  mysql(conn_id, 'close')
end

return


