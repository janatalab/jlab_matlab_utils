function outstr = mysql_is_correct(varargin)
% Checks to see whether the response to a trial is correct
%
% Input arguments are passed in as parameter/value pairs.  
% They must include:
%
% 'trial_id'
% 'session_id'
%
% By default, or if 'return_true', is set to 1, the string 'TRUE' will be
% returned if the response is correct.  If return_true=0 then 'FALSE' will be
% returned if the response is correct.  This behavior is implemented so that
% this function can service condition_matlab evaluation for contingent display
% of forms in Ensemble.

% 08/31/2009 Petr Janata

global mysql_conn_id

outstr = '';

% Parse input arguments
trial_id = [];
session_id = [];
debug_fname = '';
return_true = true;

narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'trial','trial_id'}
      trial_id = varargin{iarg+1};
    case {'session','session_id'}
      session_id = varargin{iarg+1};
    case {'debug_file'}
      debug_fname = varargin{iarg+1};
    case {'return_true'}
      return_true = varargin{iarg+1};
    otherwise
      fprintf('%s: Unknown input argument: %s\n', mfilename, varargin{iarg});
  end
end

% Check to see if we're debugging
VERBOSE = false;
if ~isempty(debug_fname)
  VERBOSE = true;
  fid = fopen(debug_fname,'wt');
  if fid < 3
    VERBOSE = false;
  end
end

if isempty(trial_id)
  error_str = sprintf('%s: Trial ID not specified', mfilename);
  if VERBOSE, fprintf(fid,error_str); end
  error(error_str);
end
if isempty(session_id)
  error_str = sprintf('%s: Session ID not specified', mfilename);
  if VERBOSE, fprintf(fid, error_str); end
  error(error_str);
end

% Check for connection to database
try 
  mysql_conn_id(1);
  conn_id = mysql_conn_id;
  if VERBOSE
    fprintf(fid, 'Using conn_id: %d\n', conn_id);
  end
  tmp_conn_id = 0;
catch   
  if VERBOSE, 
    fprintf(fid, 'No global conn_id available, making database connection\n');
  end
  mysql_make_conn;
  conn_id = 0;
  tmp_conn_id = 1;
end

% We need to get the response table name, based on the session ID
if VERBOSE, fprintf(fid,'Figuring out response table ... '); end
mysql_str = sprintf(['SELECT response_table FROM experiment WHERE ' ...
  'experiment_id = (SELECT experiment_id FROM session WHERE session_id = %d);'], session_id);
if VERBOSE, fprintf(fid,mysql_str); end
response_table = mysql(conn_id, mysql_str);
response_table = response_table{1};
if VERBOSE, fprintf(fid,'using %s\n', response_table); end

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

% Get the response to this session and trial from the response table. Make sure
% we get the most recent trial in the event that multiple iterations of this
% trial were presented in this session
mysql_str = sprintf(['SELECT response_enum, response_text, response_order FROM %s ' ...
      'WHERE session_id = %d AND trial_id = %d;'], response_table, session_id, trial_id);
[response_enum, response_text, response_order] = mysql(conn_id, mysql_str);

[dummy, resp_idx] = max(response_order);
response_text = response_text{resp_idx};

if use_enum
  if response_enum(resp_idx) == correct_response_enum
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
  
% Flip the direction of the response string if necessary
if ~return_true
  if strcmp(outstr,'TRUE')
    outstr = 'FALSE'; 
  else
    outstr = 'TRUE';
  end
end

if tmp_conn_id
  mysql(conn_id, 'close')
end

if VERBOSE, fclose(fid); end

return


