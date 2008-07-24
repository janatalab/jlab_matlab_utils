function tbl = mysql_describe_table(table, conn_id)
% Returns a description of the given mysql table.
%
% tbl = mysql_describe_table(table, conn_id);
%
% Fields are those returned by the SQL DESCRIBE command
%

% 01/26/07 Petr Janata

% Check for connection to database
try conn_id(1);
catch   
  tmp_conn_id = 1;
  mysql_make_conn;
  conn_id = 0;
end

[tbl.flds,tbl.types,tbl.null, tbl.key,tbl.default,tbl.extra] = ...
    mysql(conn_id,sprintf('DESCRIBE %s', table));

% Close the database connection if it was temporary
if (exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
