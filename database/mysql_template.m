% Stock template for setting up a connection to the mysql database on atonal
%

DEFAULT_HOST = 'atonal.ucdavis.edu';
DEFAULT_DATABASE = 'ensemble_main';

try
  host(1);
catch
  host = DEFAULT_HOST;
end;

try 
  database(1);
catch
  database = DEFAULT_DATABASE;
end

% Connect to host
mysql_make_conn(host, database);

% Close the mysql connection
mysql('close');
