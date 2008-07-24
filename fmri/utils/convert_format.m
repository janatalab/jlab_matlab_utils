function convert_format(sinfo, inpath, outpath, scan_offset) 
%  convert_format(sinfo, inpath, outpath, scan_offset);
%
%  Converts from GE format to Analyze format
%
%  sinfo is the subject info structure.  See one of the get_*_sinfo.m files for format.
%

% 06/28/00 PJ Modified to accommodate new GE_convert calling structure.
%             NOTE: Incompatible with xcult and schubert scripts
%
% 02/12/01 PJ -- implemented volume tossing during conversion step.  The
%                "starting volume" for each run will be incremented by the
%                number of volumes specified in scan_offset (default = 0) 
%
% 02/19/01 PJ -- run number assignments for functional runs increment.  This
%                prevents writing of EPI data collected in different series
%                into the same run numbers.
%
% 04/27/01 PJ Added ability to handle EPI runs distributed across multiple
%             series in a single exam. To take advantage of this specification,
%             the epi run numbers associated with a series must be specified in
%             a 3rd column of cells in the series_mapping field of the subject
%             information structure.  Otherwise, the cond_order field is used
%             and assumes that all EPI runs in an exam were collected under a
%             single series.
% 04/27/01 PJ EPI series with single and multiple runs are now handled the same way.
%             
% 04/29/01 PJ fixed indexing problem when grabbing number of volumes for a
%             run. Needed to refer to cond_order field.
% 12/04/01 PJ It turns out the s.num_vol already contains the number of volumes
%             per run ordered by run (this is set in the get_exp_sinfo.m
%             files), so there is no need to refer to cond_order.
%
% 09/04/03 PJ Added handling of double echo acquisitions in which T2 and PD
%             images are interleaved
%
% 08/30/04 PJ Allowed skipping of directories for which there are no data

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

if isempty(outpath)
  outpath = './';
end

if isempty(inpath)
  inpath = './';
end

for isub = 1:nsub
  
  s = sinfo(isub);
  subj_root = s.id;
  
  % specify root output directory and make sure it exists
  outroot = fullfile(outpath, subj_root);

  if ~exist(outroot,'dir')
    disp(sprintf('Creating directory: %s', outroot))
    unix(['mkdir ' outroot]);
  end

  num_epi = 0;
  for iexam = 1:s.nexams

    exam_root = sprintf('%05d', s.exam_nums(iexam));
    datapath = fullfile(inpath, subj_root, exam_root);
    
    series=s.series_mappings{iexam};
    nseries= size(series,1);

    for mapping_idx = 1:nseries
      indir = fullfile(datapath, char(series(mapping_idx,1)));

      % Make sure the directory whose data we want to convert exists
      if ~exist(indir)
	fprintf('Could not locate input directory: %s (skipping ...)\n', indir);
      else
	switch char(series(mapping_idx,2))
	  case {'epi','epi1','epi2','epi3','epi4','epi_12','epi_34'}
	    epidir = fullfile(outroot, char(series(mapping_idx,2)));

	    % Check to make sure that output directory exists
	    if ~exist(epidir,'dir')
	      disp(sprintf('Creating directory: %s', epidir))
	      unix(['mkdir ' epidir]);
	    end
	    
	    % Check for number of EPI runs
	    % Check to see if there is condition order information in the series
	    % mappings
	    if size(series,2) == 3
	      nruns = length(series{mapping_idx,3});
	      multi_epi_series = 1;
	    else
	      nruns = length(s.cond_order);
	      multi_epi_series = 0;
	    end
	    
	    for irun = 1:nruns
	      num_epi = num_epi+1;
	      
	      outdir = fullfile(epidir, sprintf('run%d', num_epi));

	      % Check to make sure that output directory exists
	      if ~exist(outdir,'dir')
		disp(sprintf('Creating directory: %s', outdir))
		unix(['mkdir ' outdir]);
	      end
	      
	      if multi_epi_series
		%	      nvol = s.nvol(series{mapping_idx,3}(irun));
		nvol = s.nvol(num_epi);
	      else
		nvol = s.nvol(num_epi);
	      end
	      
	      if multi_epi_series
		% 12/4/01 This code still needs to be tested before it should be
		% fully trusted
		[dummy, vol_idx] = intersect(1:length(s.cond_order), series{mapping_idx,3}(1:irun));
		run_offset = sum(s.nvol(vol_idx)) - s.nvol(vol_idx(irun));
	      else
		run_offset = sum(s.nvol(1:num_epi))-s.nvol(num_epi);
	      end

	      start_num = run_offset+2+scan_offset; % 2 = increment 1 for
	      % template, 1 to get to
	      % first image

	      outstub = fullfile(outdir, sprintf('%s_%s_r%d',subj_root, char(series(mapping_idx,1)), num_epi));
	      status = GE_convertADW(indir,outstub, start_num, ...
		  nvol-scan_offset);
	      if status
		disp(sprintf('Error in converting run %d, series %s, exam %s', irun, char(series(mapping_idx,1)), exam_root))
		% Decrement run indicator so that next run writes into the
		% location of the bad run.  Also, delete the incorrect
		% directory.
		% This functionality might be
		% assuming too much about how errors are to be handled and
		% might have to be removed in future versions
		unix(['rmdir ' outdir]);
		num_epi = num_epi - 1;
		break
	      end
	    end
	  otherwise			% e.g. {'hires','coplanar', 'coplanar2'}
	    outdir = fullfile(outroot, char(series(mapping_idx,2)));

	    outstub = fullfile(outdir, [subj_root '_' char(series(mapping_idx,1))]);
	    disp(sprintf('Source: %s; Destination: %s', indir, outdir))
	    
	    % Check to make sure that output directory exists
	    if ~exist(outdir,'dir')
	      disp(sprintf('Creating directory: %s', outdir))
	      unix(['mkdir ' outdir]);
	    end
	    switch char(series(mapping_idx,2))
	      case {'hires_T2_PD'}
		GE_convertT2PD(indir,outstub);
	      otherwise
		GE_convertVolume(indir,1,outstub);
	    end
	end % switch
      end % if exist(indir)
    end  % for mapping_idx = 1:nseries
  end % for iexam = 1:s.nexams
end % for isub = 1:nsub