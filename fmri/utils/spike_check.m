function [info, status] = spike_check(sinfo, inpath, outpath, scan_offset, nvox, sc) 
%  status = spike_check(sinfo, inpath, outpath, scan_offset, nvox, spkcrit);
%
%  Checks GE format data for spikes and dropouts in the signal.
%
%  nvox -- size of volume in voxels (x,y,z) centered on original volume to use
%  in the analysis
% 
%  spkcrit -- a structure containing fields with criteria for calling something
%  a spike:
%     .nstd - if the difference between two samples is this many standard
%             deviations (computed across the single run), the sample will be
%             flagged as bad. Default = 1.
%
%
%  The script is based on Souheil Inati's QA.m script, but modified to handle
%  multiple runs, and place various data in the output structure
%

% 04/07/01 Petr Janata
% 04/27/01 PJ Added ability to handle EPI runs distributed across multiple
%             series in a single exam. To take advantage of this specification,
%             the epi run numbers associated with a series must be specified in
%             a 3rd column of cells in the series_mapping field of the subject
%             information structure.  Otherwise, the cond_order field is used
%             and assumes that all EPI runs in an exam were collected under a
%             single series.
% 04/29/01 PJ fixed indexing problem when grabbing number of volumes for a
%             run. Needed to refer to cond_order field.
% 01/30/01 PJ Brought in line with convert_format.m (no need to check cond_order)

status = 0;

nsub = length(sinfo);

if nargin < 2
  inpath = './';
  outpath = './';
elseif nargin < 3
  outpath = './';
end

if nargin < 4
  scan_offset = 0;
end

if nargin < 5
  nvox = [15 15 2];
end

if isempty(outpath)
  outpath = './';
end

if isempty(inpath)
  inpath = './';
end

if nargin < 6
  sc.nstd = 1;
end

for isub = 1:nsub
  
  s = sinfo(isub);
  subj_root = s.id;
  
  info(isub).id = s.id;
  
  num_epi = 0;
  for iexam = 1:s.nexams

    exam_root = sprintf('%05d', s.exam_nums(iexam));
    datapath = fullfile(inpath, subj_root, exam_root);
    
    series=s.series_mappings{iexam};
    nseries= size(series,1);

    for mapping_idx = 1:nseries
      indir = fullfile(datapath, char(series(mapping_idx,1)));
      
      switch char(series(mapping_idx,2))
	case {'epi','epi1','epi2','epi3','epi4','epi_12','epi_34'}

	  % Check to see if there is condition order information in the series
	  % mappings
	  if size(series,2) == 3
	    nruns = length(series{mapping_idx,3});
	    multi_epi_series = 1;
	  else
	    nruns = length(s.cond_order);
	    multi_epi_series = 0;
	  end

	  % Create the name of the first file in inDir
	  firstfile = fullfile(indir,'I.001');

	  % Read the header
	  [su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(firstfile);

	  nX =  im_hdr.imatrix_X; % X Voxels
	  nY =  im_hdr.imatrix_Y; % Y Voxels
	  nZ =  im_hdr.slquant;   % Z Voxels
	  volSize = [nX nY nZ];

	  % Initialize
	  imageVol = zeros(nX,nY,nZ);
	  
	  xa = floor(nX/2-nvox(1));
	  xb = ceil(nX/2+nvox(1));
	  ya = floor(nY/2-nvox(2));
	  yb = ceil(nY/2+nvox(2));
	  za = floor(nZ/2-nvox(3));
	  zb = ceil(nZ/2+nvox(3));
	  nx = xb-xa+1;
	  ny = yb-ya+1;
	  nz = zb-za+1;
	  subVol = zeros(nx,ny,nz);

	  % Check for number of EPI runs
	  for irun = 1:nruns
	    num_epi = num_epi+1;
	    
	    if multi_epi_series
	      nvol = s.nvol(series{mapping_idx,3}(irun));
	    else
	      nvol = s.nvol(irun);
	    end
	    
	    info(isub).exam{num_epi} = exam_root;
	    info(isub).series{num_epi} = char(series(mapping_idx,1));
	    info(isub).run(num_epi) = num_epi;
	    
	    info(isub).data{num_epi} = zeros(nvol,1);
	    
	    if multi_epi_series
	      run_offset = sum(s.nvol(series{mapping_idx,3}(1:irun))) - s.nvol(series{mapping_idx,3}(irun));
	    else
	      run_offset = sum(s.nvol(1:irun))-s.nvol(irun);
	    end
	    
	    start_num = run_offset+1+scan_offset; % increment 1 to get to first image

	    for ivol = 1:nvol
	      [imageVol, lastfile] = GE_readVolume(indir, start_num + ivol, volSize, 16, ...
		  im_offset);
	      subVol = imageVol(xa:xb,ya:yb,za:zb);
	      info(isub).data{num_epi}(ivol) = mean(subVol(:));
	      fprintf('runnum = %d\tvolnum = %d\tread %s\n',num_epi, ivol,lastfile);
	    
	    end % for ivol=
	    
	    %
	    % Calculate some stats
	    %
	    tmp = info(isub).data{num_epi};
	    info(isub).mean(num_epi) = mean(tmp);
	    info(isub).std(num_epi) = std(tmp);
	    
	    % Find locations where there are sharp transients
	    tmpdiff = abs(diff(info(isub).data{num_epi})); 
	    badidx = find(tmpdiff > sc.nstd * info(isub).std(num_epi));
	    info(isub).badidx{num_epi} = badidx;
	    
	    info(isub).numbad(num_epi) = fix(length(badidx)/2); % **rough**
                                                                % estimate of
                                                                % number of
                                                                % spikes
								
	  
	  end % for irun=
      end % switch
    end % for mapping_idx = 
  end % for iexam = 
end % for isub=