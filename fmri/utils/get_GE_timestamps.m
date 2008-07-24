function run_info = get_GE_timestamps(sinfo,dataroot,nslice_per_vol)
% run_id_vect = get_GE_timestamps(inDir);
%
% Retrieves a vector of run ID numbers for all EPI images contained in the
% series that resides in inDir
%

% 09/27/03 PJ Added checking for presence of root directory
% 09/02/04 PJ Brought into compliance with cat syntax

nsub = length(sinfo);

if nargin < 2
  dataroot = './';
end

if nargin < 3
  error(['Please specify the number of expected slices per volume and try' ...
	' again.'])
end

for isub = 1:nsub
  s = sinfo(isub);
  subj_root = s.id;

  total_epi_series = 0;
  for iexam = 1:s.nexams
    exam_root = sprintf('%05d', s.exam_nums(iexam));
    datapath = fullfile(dataroot, subj_root, exam_root);

    series=s.series_mappings{iexam};
    epi_idx = strmatch('epi', char(series{:,2}));
    
    nser=length(epi_idx);
    if ~nser
      fprintf('Found no EPI series. Subject %s. Exam %d\n', s.id, iexam);
    else
      for iser = 1:nser
	total_epi_series = total_epi_series + 1;
	
	% Get series ID number
	inDir = fullfile(datapath, series{epi_idx(iser),1});
	ser_idstr = inDir(end);
	ser_idnum = str2num(ser_idstr);

	% Strip tailing directory
	GE_data_root = inDir(1:max(findstr('/',inDir))-1);
	
	% Make sure the root directory exists
	if ~exist(GE_data_root)
	  fprintf('Directory (%s) does not exist\n', GE_data_root);
	  break
	end

	% List directories and extract those corresponding to the desired series
	d = dir(GE_data_root);

	dir_idx = find(cat(1,d.isdir));
	ndir = length(dir_idx);

	dirmat = char(d.name);

	good_idx = find((dirmat(:,3) == ser_idstr) & isspace(dirmat(:,4)));

	ser_dirnames = char(d(good_idx).name);

	ndir = size(ser_dirnames,1);

	nimg = 0;
	im_user17 = [];
	for idir = 1:ndir
	  curdir = fullfile(GE_data_root, ser_dirnames(idir,:));
	  subd = dir(curdir);
	  
	  dirmat = char(subd.name);
	  good_idx = strmatch('I.', dirmat);
	  nfiles = length(good_idx);
	  
	  disp(sprintf('Found %d image files in directory %s', nfiles, curdir))
	  
	  flist = char(subd(good_idx).name);

	  % Loop through and read relevant field from each file.  Reading the entire
	  % header is really slow, so unless there is an explicit need to do this, don't.
	  for ifile = 1:nfiles
	    nimg = nimg+1;
	    imageFile = fullfile(curdir, flist(ifile,:));
	    
	    fid = fopen(imageFile,'r','b');
	    fseek(fid,4936,'bof');
	    im_user17(nimg) =  fread(fid,1,'float32');
	    fclose(fid);
	  end % for ifile=
	end % for idir=
	run_info.id_vect{isub,total_epi_series} = im_user17;
	run_ids = unique(im_user17);
	nrun = length(run_ids);
	for irun = 1:nrun
	  run_info.nvol{isub,total_epi_series}(irun) = sum(im_user17 == ...
	      run_ids(irun))/nslice_per_vol;
	end
      end % for iser=
    end
  end % for iexam=
end % for isub=
