function conn_id = mysql_make_conn(host, db, conn_id)
% Opens a connection to a mysql server with a specific connection ID (conn_id)
% conn_id = mysql_make_conn(host, db, conn_id);
%
% Opens a connection to a mysql server with a specific connection ID (conn_id)
% that can be used in subsequent calls to the database.  If no connection ID is
% given, a default ID of zero is used. If a connection with the specified
% conn_id is already open, it is left open.

% 01/04/07 PJ Modified to return conn_id
% 12/19/07 PJ added checking for open connection

DEFAULT_HOST = 'atonal.ucdavis.edu';
DEFAULT_DATABASE = 'ensemble_main';

try 
  host(1);
catch
  fprintf('Using default host: %s\n', DEFAULT_HOST);
  host = DEFAULT_HOST;
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

user = '';
passwd = '';
switch host
  case 	{'127.0.0.1','localhost','atonal','atonal.ucdavis.edu','atonal.cmb.ucdavis.edu'}

   switch db
    case {'ensemble_main','ensemble_tarp'}
     mysql_researcher_login; %populates user and passwd variables. This script
			     %should sit in a location that is not publically visible.
    case 'ensemble_dev'
     mysql_researcher_dev_login;
   end

end

if mysql(conn_id,'status')
  % Need to open a mysql connection
  status = mysql(conn_id,'open',host,user,passwd);

  database_str = sprintf('use %s', db);
  mysql(conn_id,database_str);
end
