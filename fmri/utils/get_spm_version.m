function version_num = get_spm_version()
% Checks to see which version of SPM5 we are dealing with. This is a bit 
% tenuous because we are assuming version.txt for SPM5 is on path before
% any other possible instance of version.txt
%
% Assume version# 958 if version.txt not found. This is really dicey.


if ~exist('version.txt')
  version_num = 957;
else
  version_info = textread('version.txt','%s');
  revidx = strmatch('Revision:', version_info);
  version_num = str2num(version_info{revidx+1});
end
