function [status,varargout] = execute_unix(unix_str,nout,errmsg,logfid)
% [varargout] = execute_unix(unix_str,nout,errmsg,logfid);
%
% Executes a command in the UNIX shell and returns results
%
%  unix_str - command to be executed
%  nout - The expected number of output arguments 1=return status only
%  errmsg - The string to be displayed if the command fails.
%  logfid - The file ID where the command should be echoed. Default is stdout.

% 03/22/06 Petr Janata

if nargin < 4
  logfid = 1;
end

if nargin < 3
  errmsg = '';
end

if nargin < 2;
  nout = 1;
end

fprintf(logfid,'%s\n', unix_str);  % Echo the command
[retval{1:nout}] = unix(unix_str);

status = retval{1};
if status
  error(errmsg)
end

for iout = 2:nout
  varargout{iout-1} = retval{iout};
end

return
