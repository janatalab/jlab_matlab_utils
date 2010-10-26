function convert_format_v2(sinfo, inpath, outpath, scan_offset) 
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
  subpath = fullfile(outpath, subj_root);
	check_dir(subpath);

	% Figure out the number of sessions we have for this subject
	nsess = length(s.sessinfo);
	
	for isess = 1:nsess
		sessinfo = s.sessinfo(isess);
		
		if ~sessinfo.use_session
			continue
		end
		
		sesspath = fullfile(subpath, sessinfo.id);
		
		% Figure out how many series we need to convert
    series = sessinfo.series_mappings;
    nseries= size(series,1);
		
		num_epi = 0;  % counter for total number of epi runs we have converted
		nexams = sessinfo.nexams;
 
		for mapping_idx = 1:nseries
			% Figure out type of series we are converting
			series_type = series{mapping_idx,2};
			
			% If we are dealing with more than one exam, figure out what the exam
			% directory is
			if nexams > 1
				exam_idx = series(mapping_idx,3);
			else
				exam_idx = 1;
			end
			
			exam_root = sprintf('%05d', sessinfo.exam_ids(exam_idx));
			exampath = fullfile(sesspath, exam_root);
    
			srcdir = fullfile(exampath, char(series(mapping_idx,1)));

      % Make sure the directory whose data we want to convert exists
      if ~exist(srcdir)
				fprintf('Could not locate input directory: %s (skipping ...)\n', srcdir);
				continue
			end
			
			switch series_type
				case {'epi','epi1','epi2','epi3','epi4','epi_12','epi_34'}
					epidir = fullfile(sesspath, series_type);
					check_dir(epidir)
					
					% Check for number of EPI runs
					% Check to see if there is condition order information in the series
					% mappings
					% In version 2 this has been moved to column 4 (epi_runcol)
					epi_runcol = 4;
					if size(series,2) == epi_runcol
						nruns = length(series{mapping_idx, epi_runcol});
						multi_epi_series = 1;
					else
						nruns = length(sessinfo.cond_order);
						multi_epi_series = 0;
					end
					
					for irun = 1:nruns
						num_epi = num_epi+1;
						
						outdir = fullfile(epidir, sprintf('run%d', num_epi));
						
						% Check to make sure that output directory exists
						check_dir(outdir)
						
						nvol = sessinfo.nvol(num_epi);
	      
						if multi_epi_series
							% 12/4/01 This code still needs to be tested before it should be
							% fully trusted
							[dummy, vol_idx] = intersect(1:length(sessinfo.cond_order), ...
								series{mapping_idx,epi_runcol}(1:irun));
							run_offset = sum(sessinfo.nvol(vol_idx)) - sessinfo.nvol(vol_idx(irun));
						else
							run_offset = sum(sessinfo.nvol(1:num_epi))-sessinfo.nvol(num_epi);
						end

						start_num = run_offset+2+scan_offset; % 2 = increment 1 for
						% template, 1 to get to first image

						outstub = fullfile(outdir, sprintf('%s_%s_r%d',sessinfo.id, char(series(mapping_idx,1)), num_epi));
						status = GE_convertADW(srcdir, outstub, start_num, nvol-scan_offset);
	      
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
					end % for irun
					
				otherwise			% e.g. {'hires','coplanar', 'coplanar2'}
					outdir = fullfile(sesspath, series_type);
					check_dir(outdir)
					
					outstub = fullfile(outdir, [subj_root '_' char(series(mapping_idx,1))]);
					disp(sprintf('Source: %s; Destination: %s', srcdir, outdir))
					
					switch series_type
						case {'hires_T2_PD'}
							GE_convertT2PD(srcdir,outstub);
						otherwise
							GE_convertVolume(srcdir,1,outstub);
					end
			end % switch % series_type
    end  % for mapping_idx = 1:nseries
  end % for isess = 1:nsess
end % for isub = 1:nsub