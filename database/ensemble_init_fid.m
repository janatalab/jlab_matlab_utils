function fid = ensemble_init_fid(params)
% Initializes a file identifier to which to direct output.
% fid = ensemble_init_fid(params);
%
% Initializes a file identifier to which to direct output.  The default in the
% absence of input parameters is standard out (fid=1).
%
% params is a structure with the following possible fields
%
% .print - a flag to indicate whether printing should occur
% .write2file - a flag to indicate whether or not printing should be done to
%               file
% .fname - name of file to print to
% .filemode - mode with which to open file. Defaults to 'wt'.

% 02/03/07 Petr Janata

% Check to see if we already have a registered output file id
if nargin > 0 && isfield(params, 'fid') && ~isempty(params.fid) && params.fid>0
  return
end

fid = 1;  % default to standard out

try print_tables = params.print; catch print_tables=1; end
try write2file = params.write2file; catch write2file = 0; end
try filemode = params.filemode; catch filemode = 'wt'; end
try verbose = params.verbose; catch verbose = 1; end

if print_tables
  
  if write2file
    fid = fopen(params.fname,filemode);
    if fid == -1
      error(sprintf('Problem opening logfile: %s\n',params.fname));
    end
    if verbose == 1
      fprintf('Writing tables to file: %s\n', params.fname);
    end
  else
    fid = 1;
  end
end
