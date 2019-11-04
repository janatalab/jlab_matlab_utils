function [conn_id, mysql_info] = mysql_make_conn(mysql_info,database,conn_id)
% Opens a connection to a mysql server.
% conn_id = mysql_make_conn(mysql_info);
%
% Opens a connection to a mysql server and returns a connection ID
% (conn_id) that can be used in subsequent calls to the database. The
% connection information is specified in the input structure (mysql_info)
% which has to have the following fields:
%
% mysql_info.host - name of host with the database
% mysql_info.database - name of database to open
% mysql_info.user - user who to connect to the database with
% mysql_info.passwd - the database password for the specific user
%
% The mysql_info structure can be obtained by calling
% mysql_login() solely with the host and database information.
%
% If no connection ID is given, a default ID of zero is used, and a conn_id
% field is attached to mysql_info and returned as an optional 2nd output argument.
% If a connection with the specified conn_id is already open, it is left open.
%
% NOTE: In contrast to earlier versions, mysql_make_conn now REQUIRES a structure
% that specifies the host, database, user, and password in order to
% establish a connection.  Failure to pass in any of this information will
% no longer result in utilization of default values.
%
% See also: mysql_login

% 01/04/07 PJ Modified to return conn_id
% 12/19/07 PJ added checking for open connection
% 10/18/09 PJ added option of passing in the first argument as a struct
%             that contains all of the information
% 10/30/09 PJ minor fix to handle empty first argument
% 11/3/09  ST minor fix (use host_or_struct in first try instead of host)
% 06/08/10 ST abstracted function to work with any host or database
%             branching based on host and database moved to mysql_researcher_login.m
% 06/15/10 PJ Modified to force passing in of connection information
  
if nargin < 1
  error('%s: Structure with following fields required:\n\t.host\n\t.database\n\t.user\n\t.passwd\n', mfilename);
end

if nargin == 3
  warning(sprintf(['The host,database,conn_id scheme for calling mysql_make_conn has been changed.\n' ...
    'The single input parameter is now a structure with the following required fields:\n' ...
    '\t.host\n\t.database\n\t.user\n\t.passwd\n\t.conn_id - for connections other than 0']))
  if ~isempty(conn_id) && ~mysql(conn_id,'status')
    fprintf('Connection with ID=%d is already open\n', conn_id)
    return
  else
    disp('Given separate login possibilities for researcher and subject, there is no way to handle your request.');
  end
end

if ~isfield(mysql_info, 'conn_id') || isempty(mysql_info.conn_id)
  mysql_info.conn_id = 0;
end

if mysql(mysql_info.conn_id,'status')
  % Need to open a mysql connection
  
  % Check to make sure we have the requisite variables to try
  if ~all(isfield(mysql_info,{'host','database','user','passwd'})) || ...
      isempty(mysql_info.host) || ...
      isempty(mysql_info.user) || ...
      isempty(mysql_info.passwd) || ...
      isempty(mysql_info.database)
    error('%s: insufficient information specified to establish connection')
  else
    conn_id = mysql(mysql_info.conn_id,'open',mysql_info.host,mysql_info.user,mysql_info.passwd);
    database_str = sprintf('use %s', mysql_info.database);
    mysql(mysql_info.conn_id,database_str);
  end
end

conn_id = mysql_info.conn_id;

return
