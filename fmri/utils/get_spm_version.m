function version_num = get_spm_version()
% Checks to see which version of SPM5 we are dealing with. This is a bit 
% tenuous because we are assuming version.txt for SPM5 is on path before
% any other possible instance of version.txt

version_info = textread('version.txt','%s');
revidx = strmatch('Revision:', version_info);
version_num = str2num(version_info{revidx+1});
