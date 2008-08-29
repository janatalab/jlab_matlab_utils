function status=run_csdm_2eeg(noffset,nsamp,avg_cells,bad_chans,ses_fname,reference)
%status=run_csdm_2eeg(noffset,nsamp,avg_cells,bad_chans,ses_fname,reference)
%version 1	1/8/96
%
%parameters
%
%	noffset = offset in samples from beginning of trial to begin analysis
%	nsamp = number of samples, must be even (converted into Npoints)
%	avg_cells = arrays of cells to use (optional)
%		supports cell collapsing, e.g. [1 2; 3 6] means 1 and 2 are
%		collapsed toproduce output cell 1 and 3 and 6 collapsed 
%		to produce output cell 2.  If no cell collapsing is wanted 
%		you must still separate output cells by semicolons
%		e.g. [ 1; 2; 3; 6];
%	bad_chans = additional bad_channels not to use (optionals)
% 	ses_fname = filename including path containing data
%	reference = reference to usein calculation (optional) defualts to 
%		the average reference 
%
%Modification History
%     Created by Ramesh Srinivasan 1/8/96

%Initialize return argument
status = -1;

%Define constants
ses_hdr_offsets_v;
TRUE = 1; FALSE = 0;

%Check number of input arguments
if ~((nargin == 6))
	error('Number of input arguments must be 6.');
end
if rem(nsamp,2) ~= 0
	error ('nsamp must be a multiple of two')
end;
if  isempty(reference)
	reference = 'average';
end; 
%Initialize input and output fids
ses_fid = -1;
ave_fid = -1;
%First try batch mode
ses_fid = fopen(ses_fname, 'r');
if ses_fid == -1
	error(['Could not open input file' ses_fname '.']);
end
%Read the session file EGIS header
[ses_fhdr,ses_chdr,ses_ename,ses_czeros,ses_cgains,ses_cnames,ses_fcom,ses_ftext,ses_coff]=rd_egis_hdr_v(ses_fid);
if isempty(avg_cells)
	avg_cells = [1:ses_fhdr(NCells)]';
end;
if isempty(noffset)
	noffset = 0;
end;
if isempty(nsamp)
	nsamp = fix(ses_chdr(1,NSamp)/2)*2;
end;
	
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

if strcmp(reference,'laplacian')
	[xelec,yelec,zelec] = electrodes(NChan+1);
end;
if strcmp(reference(1:2),'co')
	[xelec,yelec,zelec] = electrodes(NChan+1);
	sources = xyz2tp(xelec,yelec,zelec);
	A = transfer_matrix(500,0.2,80,8,8.2,8.7,9.2,7.8,sources,sources);
end;
num_cells_avg = size(avg_cells,1);
%Begin csdm calculation	
for c=1:num_cells_avg
	accessed_cell = 0;
	ncol = size(find(avg_cells(c,:)),2);
	disp(['Processing output cell ' int2str(c)]);
	accessed_cell = FALSE;
	start_samp = noffset +1;
	stop_samp = noffset +nsamp;
	ncsdm = (ses_fhdr(NChan)+2)*(ses_fhdr(NChan)+1)/2;
	avgdata = zeros(nsamp/2,ncsdm);
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
				trialdata2 = cal_gain(trialdata,ses_cgains,ses_czeros);
				%Reference Matrix
				if ~(strcmp(reference,'laplacian')|strcmp(reference(1:2),'co'))
					[ref_trialdata,masknew] = new_reference(trialdata2,ses_mask(ses_mask_cell_offset + t,:),reference);
					ref_trialdata = ref_trialdata.*(ones(size(ref_trialdata,1),1)*masknew);
					ses_mask(ses_mask_cell_offset+t,:) = masknew;
				elseif strcmp(reference,'laplacian')
					good_chan = find(ses_mask(ses_mask_cell_offset + t,:));	
					ref_trialdata = laplacian_trial([trialdata2 zeros(size(trialdata2,1),1)],good_chan,xelec,yelec,zelec);
					masknew = ses_mask(ses_mask_cell_offset + t,:);
				elseif strcmp(reference(1:8),'cortical')
					good_chan = find(ses_mask(ses_mask_cell_offset + t,:));
					good_trial = zeros(size(trialdata2,1),size(good_chan,2));
					trialdata3 = average_reference(trialdata2(:,1:NChan),ses_mask(ses_mask_cell_offset + t,:));
					good_trial = trialdata3(:,good_chan);
					good_A = zeros(size(good_chan,2),size(good_chan,2));
					good_A = A(good_chan,good_chan);
					ref_trialdata = zeros(size(trialdata3,1),size(trialdata3,2));
					sigma_m = 1;
					sigma_v = str2num(reference(9:10));
			
					[m,deviations] = bayes_dipole_trial(good_A,good_trial',sigma_v,sigma_m);
					ref_trialdata(:,good_chan) = m';
					masknew = ses_mask(ses_mask_cell_offset + t,:);
				end
					
				%zero average epoch
				ref_trialdata = zeromean(ref_trialdata,start_samp,stop_samp);					
				%add to average stack
					
				avgdata = avgdata + csdm(ref_trialdata(start_samp:stop_samp,:));
				%else
				%disp(['Skipping trial ' int2str(t)]);
			end
		end
	end;
	num_good_trials = zeros(1,ses_fhdr(NChan)+1);
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
	avgdata = avgdata./(ones(nsamp/2,1)*good_cross);
	disp(['Writing csdm data file output cell:' int2str(c)]);
	ave_fname = [ses_fname(5:12) 'tc' int2str(c) '.' reference(1:4) '.csdm'];
	ave_fid = fopen(ave_fname,'wb'); 
	version = -1;
	fwrite(ave_fid,version,'int16');
	SampRate
	ses_chdr(c,SampRate)
	Epoch = nsamp/ses_chdr(c,SampRate);
	fwrite(ave_fid,Epoch,'int16');
	fwrite(ave_fid,ses_chdr(c,NTrials),'int16')
	fwrite(ave_fid,max(good_cross),'int16');
	Nbad_chan = size(bad_chans,2);
	fwrite(ave_fid,Nbad_chan,'int16');
	fwrite(ave_fid,bad_chans,'int16');
	if strcmp(reference,'average')
		ref_flag = 1;
	elseif strcmp(reference,'avgmast')
		ref_flag = 2;
	elseif strcmp(reference,'perimeter');
		ref_flag = 3;
	elseif strcmp(reference,'vertex');
		ref_flag = 4;
	elseif strcmp(reference,'laplacian');
		ref_flag = 5;
	elseif strcmp(reference(1:8),'cortical');
		ref_flag = 7;
	else
		ref_flag = 6;
	end;
	fwrite(ave_fid,ref_flag,'int16');
	if ref_flag == 6
		fwrite(ave_fid,size(reference,2),'int16');
		fwrite(ave_fid,reference,'int16');
	end
	if ref_flag == 7
		sigma_v = str2num(reference(9:10));
		fwrite(ave_fid,sigma_v,'int16');
	end;
	max_freq = 50*Epoch + 1;
	fwrite(ave_fid,max_freq,'int16');
	fwrite(ave_fid,size(avgdata,2),'int16');
	num_written = fwrite(ave_fid, (real(avgdata(1:max_freq,:)))', 'float');
	num_written = fwrite(ave_fid, (imag(avgdata(1:max_freq,:)))', 'float');
	fclose(ave_fid);
end

fclose('all');
disp('Finished CSDM Calculation.');
status = 1;






