function [fld_vals] = mysql_unique_fld_vals(tbl_name, fld_name, conn_id)
% Returns a sorted list of unique values in a specific field of a specific table
%
% [fld_vals] = mysql_unique_fld_vals(tbl_name, fld_name, conn_id);
%

% 08/18/05 Petr Janata
% 01/03/07 PJ - removed dependency on unique()

% Do some input parameter checking. 
% Modify min_arg and max_arg accordingly for your function
min_arg = 2;
max_arg = 3;

msg = nargchk(min_arg,max_arg,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

% Connect to host with a temporary connection if necessary
try conn_id(1);
catch   
  tmp_conn_id = 1;
  mysql_make_conn;
  conn_id = 0;
end

mysql_str = sprintf('SELECT DISTINCT %s FROM %s', fld_name, tbl_name);
fld_vals = mysql(conn_id, mysql_str);

%
% Close the mysql connection if this was a temporary opening of the database
%

if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
