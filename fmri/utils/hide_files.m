function [] = hide_files(fi)
% hide_files(fi)
%
% fi.rootdir = root directory that contains subject directories
% fi.destdir = name stub of directory to move files into
% fi.prestr = cell array of file prefixes to process
% fi.sub_id{sub1, sub2} = cell array of subject ids
% fi.vol_info(1:nsub) = array of vol_info structures
%    vol_info.run_id() = array of runs numbers to use
%    vol_info.vol_id{1:nrun} = vectors of volumes for each run to move to destdir
%

nsub = length(fi.sub_id);

for isub = 1:nsub
  disp(sprintf('Moving files for subject %s', fi.sub_id{isub}))
  subpath = fullfile(fi.rootdir,fi.sub_id{isub});
  
  for irun = 1:length(fi.vol_info(isub).run_id)
    run_idx = fi.vol_info(isub).run_id(irun);
    disp(sprintf('    Run %d', run_idx))
    
    srcdir = fullfile(subpath,sprintf('epi/run%d',run_idx));
    destdir = fullfile(subpath,sprintf('epi/run%d',run_idx),fi.destdir);
    
    if ~exist(destdir,'dir')
      disp(sprintf('Creating directory: %s', destdir))
      unix(['mkdir ' destdir]);
    end

    nvol = length(fi.vol_info(isub).vol_id{irun});
    for ivol = 1:nvol
      vol_idx = fi.vol_info(isub).vol_id{irun}(ivol);
      
      for ipre = 1:length(fi.prestr)
	unix_str = sprintf('mv %s/%s%s*i%04d.* %s', srcdir, fi.prestr{ipre},fi.sub_id{isub}, vol_idx, destdir);
	unix(unix_str);
      end
    end % for ivol =
  end % for irun=
end % for isub=
