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

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

% INSERT FUNCTION SPECIFICS HERE

return
