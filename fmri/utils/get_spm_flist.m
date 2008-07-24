function flist = get_spm_flist(srcdir,srcstub)
% flist = get_spm_flist(srcdir);
%
% Gets listing of .img files in a directory and returns a cell array suitable
% for passing to SPM5 routines

% 02/14/06 Petr Janata
  
if nargin < 2
  srcstub = '*.img';
end
  
dirlist = dir(fullfile(srcdir,srcstub));
nfiles = length(dirlist);
flist = [repmat(srcdir,nfiles,1) repmat('/',nfiles,1) ...
      strvcat(dirlist.name)];  %  repmat(',1',nfiles,1)
%flist = cellstr(flist);