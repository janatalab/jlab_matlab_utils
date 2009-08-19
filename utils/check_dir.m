function status = check_dir(outdir, verbose, parents)
% Checks for existence of outdir, and creates it if necessary
%
% status = check_dir(outdir, verbose, recurse);
%
% status returns 0 on success and 1 on failure
%
% INPUT
% outdir:  a string that specifies the directory to check.
% verbose: 0 or 1. 1 means that some basic status messages will
%          output to the command line
% parents: 0 or 1. 1 means that check_dir will execute mkdir with the '-p'
%          argument, if outdir does not exist. This will create parent
%          directories as needed. execute 'man mkdir' from the unix command
%          line for more information
%
% Copyright (c) 2005 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% 10/18/05 Petr Janata
% 08/17/09 Fred Barrett - added recursive option

if ~exist('verbose','var')
  verbose = 0; 
end
if exist('parents','var') && parents
  mkdiropts = '-p ';
else
  mkdiropts = '';  
end
status = -1;

if exist(outdir) ~= 7
  if verbose
    fprintf('Making directory: %s\n', outdir);
  end
  %escape spaces and other strange characters for mkdir
  outdir = regexprep(outdir,'[()'',& ]','\\$0');
    
  unix_str = sprintf('mkdir %s%s',mkdiropts,outdir);
  status = unix(unix_str);

else
  status = 0;
end
