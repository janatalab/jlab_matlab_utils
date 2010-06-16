function sessinfo = mysql_get_sessinfo(varargin)
% Returns session information for given sessions.
% sessinfo = mysql_get_sessinfo(varargin)
%
% Returns information associated with one or more sessions provided in an array
% of sessions. Arguments are passed in as tag/value pairs.
%
% Supported input arguments (tags):
% 'session_id' - vector of session IDs
% 'conn_id' - mysql connection ID to utilize - REQUIRED

% 01/26/07 Petr Janata
% 06/01/07 Stefan Tomic - fixed handling of temp conn_ids
% 06/10/15 PJ mysql_make_conn

% Initialize some variables
sessinfo = ensemble_init_data_struct;

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'session_id','session_ids'}
      sessids = varargin{iarg+1};
    case 'conn_id'
      conn_id = varargin{iarg+1};
    otherwise
      fprintf('mysql_get_sessinfo: Unknown input argument: %s\n', varargin{iarg});
  end
end

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

sessinfo.type = 'session_info';

% Call mysql_extract_data
[sessinfo.data, sessinfo.vars] = mysql_extract_data('table','session','session_id',sessids,'conn_id',conn_id);

return