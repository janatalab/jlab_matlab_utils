function status = check_dir(outdir, verbose)
% Checks for existence of outdir, and creates it if necessary
%
% status = check_dir(outdir, verbose);
%
% status returns 0 on success and 1 on failure
%
% INPUT
% outdir:  a string that specifies the directory to check.
% verbose: 0 or 1. 1 means that some basic status messages will
%          output to the command line
%
% Copyright (c) 2005 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% 10/18/05 Petr Janata

if ~exist('verbose','var')
  verbose = 0; 
end
status = -1;

if exist(outdir) ~= 7
  if verbose
    fprintf('Making directory: %s\n', outdir);
  end
  %escape spaces and other strange characters for mkdir
  outdir = regexprep(outdir,'[()'',& ]','\\$0');
    
  unix_str = sprintf('mkdir %s', outdir);
  status = unix(unix_str);

else
  status = 0;
end
