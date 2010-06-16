function tbl = mysql_describe_table(table, conn_id, params)
% Returns a description of the given mysql table.
%
% tbl = mysql_describe_table(table, conn_id);
%
% Fields are those returned by the SQL DESCRIBE command
%

% 01/26/07 Petr Janata
% 06/15/10 PJ sanitized mysql_make_conn

% Check for connection to database
try conn_id(1);
catch   
  if exist('params','var')
    tmp_conn_id = 1;
    if ~isfield(params,'conn_id')
      params.conn_id = 0;
    end
    conn_id = mysql_make_conn(params);
  else
    error('%s: No valid conn_id provided', mfilename)
  end
end

[tbl.flds,tbl.types,tbl.null, tbl.key,tbl.default,tbl.extra] = ...
    mysql(conn_id,sprintf('DESCRIBE %s', table));

% Close the database connection if it was temporary
if (exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
