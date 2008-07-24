function new_stimlist = check_stim_dirs(stimlist,varargin)
% Creates a directory structure for ipem_proc_series calculations 
%
% new_stimlist = check_stim_dirs(stimlist);
%
% Given a list of paths to stimuli, checks to make sure that there are
% corresponding directories in the destination tree.  At the level of the
% specific stimulus file, a directory is created with the stimulus name.  A
% subdirectory is placed in that directory with the extension for the
% filetype. Within that directory, a symlink is optionally created to the
% original stimulus file. 
%
% Returns a list of files with the new target path.
%
% Subsequent analysis directories will be placed at this same level.
%
% INPUTS
% stimlist - relative paths from srcroot to individual stimuli, passed as a cell array of strings
% srcroot  - (optional) passed in as tag/value pair. Directory where stimuli are stored.
% destroot - (optional) passed in as tag/value pair. Directory where subdirectories for each stim will be created
%
% OUTPUTS
% new_stimlist - list of symlinks to original stimuli.
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved
%
% December 29, 2006 Petr Janata
% April 30,2008 - Added escaping of irregular directory name
%                 characters(space, parmenthesis, etc.)

destroot = '';
srcroot = '';
MAKE_SYMLINK = false;
VERBOSE = false;
new_stimlist = {};

% Process the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'destroot'}
      destroot = varargin{iarg+1};
    case {'srcroot'}
      srcroot = varargin{iarg+1};
      MAKE_SYMLINK = true;
    case {'verbose'}
      VERBOSE = varargin{iarg+1};
    otherwise
      fprintf('check_stim_dirs: Unknown input argument: %s\n', varargin{iarg});
  end
end

if ~iscell(stimlist)
  stimlist = {stimlist};
end


if(VERBOSE)
  fprintf('\nChecking stimulus directories...\n');
end


nstims = length(stimlist);
for istim = 1:nstims
  if VERBOSE
    fprintf('Stim %d/%d\n', istim,nstims);
  end
  curr_dir = destroot;
  res = stimlist{istim};
  while ~isempty(res)
    [tok,res] = strtok(res,filesep);
    if ~isempty(res)
      curr_dir = fullfile(curr_dir,tok);
      check_dir(curr_dir, VERBOSE);
    else
      [dummy,fname,fext] = fileparts(tok);
      curr_dir = fullfile(curr_dir,fname);
      check_dir(curr_dir, VERBOSE);  % directory named after the stimulus
      
      curr_dir = fullfile(curr_dir,fext(2:end));
      check_dir(curr_dir, VERBOSE);
      
      if MAKE_SYMLINK
	srcfile = fullfile(srcroot,stimlist{istim});
	linkname = fullfile(curr_dir,tok);
	if ~exist(linkname)
      srcfile = regexprep(srcfile,'[()'',& ]','\\$0');
      linkname_escaped = regexprep(linkname,'[()'',& ]','\\$0');
	  unix_str = sprintf('ln -s %s %s >& /dev/null', srcfile, linkname_escaped);
	  if VERBOSE
	    fprintf('%s\n', unix_str);
	  end
	  status = unix(unix_str);
      
	end
	new_stimlist(end+1) = {linkname};
      end
    end
  end
end % for istim=
