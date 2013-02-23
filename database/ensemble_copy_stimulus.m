function status = ensemble_copy_stimulus(stiminfo,params)
% Copies a stimulus from the location specified by the combination of 
% params.stimroot and params.location to the destination directory in
% params.outpath
%
% Intended for use in conjunction with ensemble_processStims

% 23Feb2013 Petr Janata


status = 0;

srcfname = fullfile(params.stimroot,params.location);
if ~exist(srcfname,'file')
  fprintf('Could not locate source stimulus: %s\n', srcfname);
  return
end

[fpath,fname,fext] = fileparts(params.location);
destfname = fullfile(params.outpath, [fname fext]);

unix_str = sprintf('cp %s %s', srcfname, destfname);
status = ~unix(unix_str);

return
end