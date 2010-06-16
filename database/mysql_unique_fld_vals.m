function [fld_vals] = mysql_unique_fld_vals(tbl_name, fld_name, conn_id)
% Returns a sorted list of unique values in a specific field of a specific table
% [fld_vals] = mysql_unique_fld_vals(tbl_name, fld_name, conn_id);
%
% conn_id - connection to database - required

% 08/18/05 Petr Janata
% 01/03/07 PJ - removed dependency on unique()
% 06/15/10 PJ - sanitized mysql_make_conn()


% Do some input parameter checking. 
% Modify min_arg and max_arg accordingly for your function
min_arg = 2;
max_arg = 3;

msg = nargchk(min_arg,max_arg,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

mysql_str = sprintf('SELECT DISTINCT %s FROM %s', fld_name, tbl_name);
fld_vals = mysql(conn_id, mysql_str);
