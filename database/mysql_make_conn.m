function conn_id = mysql_make_conn(host_or_struct, db, conn_id)
% Opens a connection to a mysql server with a specific connection ID (conn_id)
% conn_id = mysql_make_conn(host_or_struct, db, conn_id);
%
% Opens a connection to a mysql server with a specific connection ID (conn_id)
% that can be used in subsequent calls to the database.  If no connection ID is
% given, a default ID of zero is used. If a connection with the specified
% conn_id is already open, it is left open.
%
% The first input argument can be a structure that contains the host,
% database, and conn_id parameters.
%

% 01/04/07 PJ Modified to return conn_id
% 12/19/07 PJ added checking for open connection
% 10/18/09 PJ added option of passing in the first argument as a struct
%             that contains all of the information
% 10/30/09 PJ minor fix to handle empty first argument
% 11/3/09  ST minor fix (use host_or_struct in first try instead of host)
% 06/08/10 ST abstracted function to work with any host or database
%             branching based on host and database moved to mysql_researcher_login.m
  
%ensemble_globals loads DEFAULT_HOST and DEFAULT_DATABASE
%see ensemble_globals_template.m for further information
ensemble_globals;

try 
  host_or_struct(1);
catch
  fprintf('Using default host: %s\n', DEFAULT_HOST);
  host_or_struct = DEFAULT_HOST;
end

if isstruct(host_or_struct)
  try 
    host = host_or_struct.host;
  catch
    host = DEFAULT_HOST; 
  end
  
  try
    db = host_or_struct.database;
  catch
    db = DEFAULT_DATABASE;
  end
  
  try
    conn_id = host_or_struct.conn_id;
  catch
    conn_id = 0;
  end
elseif isempty(host_or_struct)
  host = DEFAULT_HOST;
else
  host = host_or_struct;
end

try 
  db(1);
catch
  fprintf('Using default database: %s\n', DEFAULT_DATABASE);
  db = DEFAULT_DATABASE;
end

try conn_id(1);
catch conn_id = 0;
end;

%populates user and passwd
%see mysql_researcher_login_template.m for further information
ensemble = struct('host',host,'database',db);
ensemble = mysql_researcher_login(ensemble);

if mysql(conn_id,'status')
  % Need to open a mysql connection
  status = mysql(conn_id,'open',ensemble.host,ensemble.user,ensemble.passwd);

  database_str = sprintf('use %s', db);
  mysql(conn_id,database_str);
end
