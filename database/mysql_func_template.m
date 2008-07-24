% MySQL connection and argument checking utility.
% mysql_func_template.m
%
% function outargs = myfunc(varargin)
%
% Template script for checking input arguments for a connection ID and
% establishing a connection if necessary.
%
% This is good for various utility scripts that are used to retrieve specific
% pieces of information from the database
%

% Do some input parameter checking. 
% Modify min_arg and max_arg accordingly for your function
min_arg = 1;
max_arg = 2;

msg = nargchk(min_arg,max_arg,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case 'conn_id'
      conn_id = varargin{iarg+1};
    otherwise
      fprintf('get_qtxt: Unknown input argument: %s\n', varargin{iarg});
  end
end

% Check for connection to database
try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
end

% INSERT FUNCTION SPECIFICS HERE

if ~conn_id
  mysql(conn_id,'close');
end

