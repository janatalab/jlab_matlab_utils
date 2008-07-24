function clean_exp_files(sinfo, modtype_idx, modlist_idx, rootdir)
% clean_exp_mfiles.m
%
% Performs modifications on files in epi directories
%
% Specify list of indices for things to do (modtype)
% and types of files to modify (modlist)
%
% modtype = {'chmod 444', 'rm'}
% modlist = {'GE','orig_adw','normed','snormed','realigned','srealigned','mean'}
%

% 02/23/01 PJ

error(nargchk(3,4,nargin))

if nargin < 4
  rootdir = './';
end

modtype = {'chmod 444', 'rm', 'ls'};

modlist = {'GE','orig_adw','normed','snormed','realigned','srealigned','mean'};

% Get subject info
nsub = length(sinfo);

modstr = modtype{modtype_idx};
nmod = length(modlist_idx);
for isub = 1:length(sinfo)
  for imod = 1:nmod
    fpath = fullfile(rootdir, sinfo(isub).id, 'epi');
    
    switch modlist{modlist_idx(imod)}
      case 'orig_adw'
	prefix = sinfo(isub).id;
      case 'normed'
	prefix = 'n';
      case 'snormed'
	prefix = 'sn';
      case 'realigned'
	prefix = 'r';
      case 'srealigned'
	prefix = 'sr';
      case 'mean'
	prefix = 'mean';
      otherwise
	disp(sprintf('Do not know how to modify: %s', modlist(modlist_idx(imod))))
    end

    % Delete run by run to circumvent list length limits
    d = dir(fpath);
    dirlist = strvcat(d(strmatch('run',strvcat(d.name))).name);
    for idir = 1:size(dirlist,1)
      unix_str = sprintf('%s %s/run%d/%s*', modstr, fpath, idir, prefix);
      unix(unix_str);
    end
  end
end