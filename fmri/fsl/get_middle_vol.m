function volidx = get_middle_vol(fname, logfid)
% volidx = get_middle_vol(fname, logfid);
%
% Returns the index of the middle volume in a 4D data volume
%

% 10/26/05 Petr Janata

% check existence of file
if ~exist(fname,'file') && ~exist([fname '.nii'],'file') ...
        && ~exist([fname '.nii.gz'],'file')
  error('get_middle_vol:Could not find file: %s', fname)
end

try logfid;
catch logfid = 1;
end

fsl_str = sprintf('fslnvols %s', fname);
fprintf(logfid,'%s\n', fsl_str);
[status, volidx] = unix(fsl_str);  % get the number of volumes
volidx = str2num(volidx);
volidx = fix(volidx/2);  % take middle timepoint

if volidx == 0
  volidx = 1;
end
