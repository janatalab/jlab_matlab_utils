function status = COG_startup(bel_matlab_path)
%status = COG_startup(bel_matlab_path)
%
%gives you access to standard m files available for PPC on COG
%
%written by RS on 10/9/95
%
%modified by PJ on 1/18/96 -- user may now specify path leading up to bel_matlab
%structure
%		Script should now be portable accross the machines we use without
%needing to
%		modify the entire path structure each time
%

%  check to see if user has specified a path
exist('bel_matlab_path')
if ~exist('bel_matlab_path')
  if findstr(computer,'MAC')
    bel_matlab_path = 'COG:Shared:MATLAB:';		% default to COG directory if
%Mac is detected
  else
    bel_matlab_path = '/net/hebb/external2/ext2/users/exports/matlab/';	%
%assume user is on Darkwing
  end
end

path(path, bel_matlab_path);
path(path,[bel_matlab_path 'egis_file_utils']);
path(path,[bel_matlab_path 'data_utils']);
path(path,[bel_matlab_path 'scripts']);
path(path,[bel_matlab_path 'misc_utils']);
path(path,[bel_matlab_path 'ppc_utils']);
path(path,[bel_matlab_path 'net_utils']);
path(path,[bel_matlab_path 'meg_code']);
path(path,[bel_matlab_path 'spline']);
path(path,[bel_matlab_path 'phantom']);
path(path,[bel_matlab_path 'permute']);
status = 1;







