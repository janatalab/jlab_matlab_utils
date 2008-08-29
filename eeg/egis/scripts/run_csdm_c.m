function status=run_csdm_c(noffset,nsamp,nwin,avg_cells,bad_chans,ses_fname,ave_fname)
%status=run_csdm(noffset,nsamp,nwin,avg_cells,bad_chans,ses_fname,ave_fname)%
%version 1	1/8/96
%
%parameters
%
%	noffset = offset in samples from beginning of trial to begin analysis
%	nsamp = number of samples, must be even (converted into Npoints)
%   	nwin = # of windows (converted into NObs)run_csdm_c
%	avg_cells = arrays of cells to use (optional)
%	bad_chans = additional bad_channels not to use (optionals)
% 	ses_fname = filename including path containing data
%	ave_fname = filename containing csdm matrices. 
%
%Modification History
%     Created by Ramesh Srinivasan 1/8/96

%Initialize return argument
status = -1;

%Define constants
ses_hdr_offsets_v;
ave_hdr_offsets_v;
TRUE = 1; FALSE = 0;

%Check number of input arguments
if ~((nargin == 3)|(nargin == 5)|(nargin == 7))
	error('Number of input arguments must be either 3, 5, or 7.');
end
if rem(nsamp,2) ~= 0
	error ('nsamp must be a power of two')
end;
%Initialize input and output fids
ses_fid = -1;
ave_fid = -1;
%First try batch mode
if nargin == 7
	ses_fid = fopen(ses_fname, 'r');
	if ses_fid == -1
		error(['Could not open input file' ses_fname '.']);
	end
	ave_fid = fopen(ave_fname, 'wb');
	if ave_fid == -1
		temp =fclose(ses_fid);
		error(['Could not open output file ' ave_fname '.']);
	end
%Otherwise run interactive mode
else
	while ses_fid == -1
		[ses_fid, ses_fname]=get_fid('r','*.ses*', 'Open Session File:');
	end
	while ave_fid == -1
		[ave_fid, ave_fname]=put_fid('wb','.csdm','Save New csdm File As:');
	end
end

%Read the session file EGIS header
[ses_fhdr,ses_chdr,ses_ename,ses_czeros,ses_cgains,ses_cnames,ses_fcom,ses_ftext,ses_coff]=rd_egis_hdr_v(ses_fid);

%Assign defaults based on header info
if nargin == 3
	bad_chans = [];
	avg_cells = [1:ses_fhdr(NCells)];
end
if ((nargin == 5)&(length(avg_cells) == 0))
	avg_cells = [1:ses_fhdr(NCells)];
end

%More error checking. Compare arguments to header info.
%Too many cell?
%if length(avg_cells) > ses_fhdr(NCells)
% 	error(['Too many cells specified for this data file.']);
%end
%Bullshit cell number?
if any(avg_cells > ses_fhdr(NCells))
	error('A specified cell number is greater than the number of cells in file.');
end
%Too many bad_chans?
if length(bad_chans) > ses_fhdr(NChan)
 	error(['Too many bad channels specified for this data file.']);
end
%Check bad_chans argument
if any(bad_chans > ses_fhdr(NChan))
	error('A specified bad channel is greater than the number of channels in file.');
end



%Read the MacAverager .edc edit codes file
ses_mask = artifact_edit(ses_fname, ses_fhdr, ses_chdr,0);
%add the bad channels mask argument to the .edc information
ses_mask(:, bad_chans) = zeros(size(ses_mask,1),length(bad_chans));

%Construct the new output average file header
disp('Constructing the average header');
% Copy the session file header
num_cells_avg = size(avg_cells,1);
ave_fhdr = ses_fhdr; 
ave_fhdr(NChan) = (ses_fhdr(NChan)+1)*(ses_fhdr(NChan) + 2)/2;
ave_chdr(1:num_cells_avg,:) = ses_chdr(avg_cells(:,1),:); 
ave_ename = ses_ename; ave_fcom = ses_fcom; ave_ftext = ses_ftext;
ave_fhdr(NCells) = num_cells_avg;
ave_fhdr(LastDone) = 1;
ave_cnames(1:num_cells_avg,:) = ses_cnames(avg_cells(:,1),:);

% Copy the experiment info from the session file to the average file tspecs
[ave_chdr] = sesfhdr2avetspecs(ses_fhdr, ave_chdr);

%Determine the number of good trials in the average
for c = 1:num_cells_avg
	ave_chdr(c,NAvg) = 0;		% Zero the number of good trials
	ave_chdr(c,NObs) = nwin;
	ave_chdr(c,NPoints) = nsamp;
	cell_offset_mask(c) = sum(ses_chdr(1:avg_cells(c,1), NTrials)) - ses_chdr(avg_cells(c,1), NTrials);
	for t=1:ses_chdr(avg_cells(c,1),NTrials)
		trial_offset = cell_offset_mask(c) + t;
		if sum(ses_mask(trial_offset,:)) > 0
			ave_chdr(c, NAvg) = ave_chdr(c, NAvg) + 1;
		end
	end
end


% Write out the average file header
wt_csdm_hdr_v(ave_fid,ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext);
fseek(ave_fid, 0, 'eof');



%Begin csdm calculation	
for c=1:num_cells_avg
	accessed_cell = 0;
	ncol = size(find(avg_cells(c,:)),2)
	disp(['Processing output cell ' int2str(c)]);
	for iwin = 1:nwin
		disp(['  Observation ' int2str(iwin)]);
		accessed_cell = FALSE;
		start_samp = noffset +1 +(iwin-1)*nsamp;
		stop_samp = noffset +iwin*nsamp;
		avgdata = zeros(ave_chdr(c,NPoints)/2, ave_fhdr(NChan));
		power_trialdata = zeros(ave_chdr(c,NPoints)/2, ave_fhdr(NChan));
		for icol = 1:ncol
		thecell = avg_cells(c,icol);
			for t=1:ses_chdr(thecell, NTrials)
				ses_mask_cell_offset = sum(ses_chdr(1:thecell, NTrials)) - ses_chdr(thecell,NTrials);
				if sum(ses_mask(ses_mask_cell_offset + t,:)) > 0	
					if accessed_cell == FALSE
						trialdata = rd_onetr_allch(ses_fid, ses_coff(thecell), t, ses_fhdr(NChan), ses_chdr(thecell, NPoints), 'bof');
						accessed_cell = TRUE;
					else
						trialdata = rd_onetr_allch(ses_fid, ses_coff(thecell), t, ses_fhdr(NChan), ses_chdr(thecell, NPoints), 'cof');
					end
				%Convert matrix to microvolts
					vol_trialdata = cal_gain(trialdata,ses_cgains,ses_czeros);
				%Average Reference Matrix
					avref_trialdata= average_reference(vol_trialdata, ses_mask(ses_mask_cell_offset + t,:));
				%zero average epoch
					avref_trialdata = zeromean(avref_trialdata,start_samp,stop_samp);			
				%	call to csdm
					power_trialdata = csdm(avref_trialdata(start_samp:stop_samp,:));		
				%add to average stack
					avgdata = avgdata + power_trialdata;
				%else
				%disp(['Skipping trial ' int2str(t)]);
				end
			end
		end;
		num_good_trials = zeros(1,ses_fhdr(NChan)+1)
		for icol = 1:ncol
			thecell = avg_cells(c,icol);
			ses_mask_cell_offset = sum(ses_chdr(1:thecell,NTrials)) - ses_chdr(thecell,NTrials);
			start_offset = ses_mask_cell_offset + 1;
			stop_offset = start_offset + ses_chdr(thecell,NTrials) - 1;
			num_good_trials = num_good_trials+ sum(ses_mask(start_offset:stop_offset,:));
		end;
		%divides are different for amplitude and phase
		icount = 1;
		good_cross = zeros(1,size(avgdata,2));
		for ichan = 1:ses_fhdr(NChan)+1
			for jchan = ichan:ses_fhdr(NChan)+1
				good_cross(icount) = min([num_good_trials(ichan) num_good_trials(jchan)]);
				if good_cross(icount) == 0
					good_cross(icount) = 1;
				end;
				icount = icount+1;
			end;
		end;
		avgdata = avgdata./(ones(ave_chdr(c,NPoints)/2,1)*good_cross);
		disp('Writing data');
		num_written = fwrite(ave_fid, (real(avgdata))', 'float');
		num_written = fwrite(ave_fid, (imag(avgdata))', 'float');
		if thecell==1 & iwin == 2
			save avgcsdm.mat avgdata
		end;
		if num_written ~= size(avgdata,1) * size(avgdata,2)
			error(['Failed to write data of cell ' int2str(c)]);
		end
	end;
end

fclose('all');
disp('Finished CSDM Calculation.');
status = 1;


