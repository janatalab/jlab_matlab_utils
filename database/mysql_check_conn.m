function status = mysql_check_conn(conn_id,action)
% Checks a mysql connection, opens or closes if necessary.
% conn_id = mysql_check_conn(conn_id,action);
%
% Checks to see if a mysql connection with the desired ID is active and checks
% to see if there is an action that should be performed.

% 06/15/10 PJ Sanitized mysql_make_conn

if nargin < 2
  action = 'check';
end

status = mysql(conn_id,'status'); % returns 0 if active

if status==0 && ~strcmp(action,'close')
  return
end

switch action
  case 'close'
    mysql(conn_id,'close');
    status = 1;
end
