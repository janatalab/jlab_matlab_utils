function sessinfo = mysql_get_sessinfo(varargin)
% Returns session information for given sessions.
% sessinfo = mysql_get_sessinfo(varargin)
%
% Returns information associated with one or more sessions provided in an array
% of sessions. Arguments are passed in as tag/value pairs.
%
% Supported input arguments (tags):
% 'session_id' - vector of session IDs
% 'conn_id' - mysql connection ID to utilize

% 01/26/07 Petr Janata
% 06/01/07 Stefan Tomic - fixed handling of temp conn_ids

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

% Check for connection to database
try conn_id(1);
catch   
  tmp_conn_id = 1;
  mysql_make_conn;
  conn_id = 0;
end

sessinfo.type = 'session_info';

% Call mysql_extract_data
[sessinfo.data, sessinfo.vars] = mysql_extract_data('table','session','session_id',sessids,'conn_id',conn_id);

% Close the database connection if it was temporary
if (exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
