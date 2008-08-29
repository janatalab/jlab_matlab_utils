function status = run_meg_csdm(skip_samples,Epoch,MEG_max,EEG_max,EEG_reference,add_bad_chan,max_freq,megfilename);
%status = run_meg_csdm(skip_samples,Epoch,MEG_max,EEG_max,EEG_reference,
%	add_bad_chan,max_freq,meg_filename);
%
%skip_samples = # of samples from beginning of file to average 
%Epoch = epoch length in seconds to use
%MEG_max = maximum value (in picotesla) to keep an MEG channel
%EEG_max = maximum value (in microvolts) to keep an EEG channel
%EEG_reference = reference to use for EEG (defaults to average mastoids)
%add_bad_chan = bad_channels to ignore
%max_freq = highest frequency to write out 
%megfilename = meg_filename (must end in .MEG or .meg) gui o.k. by skipping
% channels are sequenced in the following absolute channel order:
% 1 to 148 :: MEG
% 149 to 149 + NReference :: reference channels  
% 149 +NReference +1 to 149 +NReference + NEeg :: EEG signals 
%	Note:  this is not necessarily the same in the csdm file 
%		it depends on the EEG reference selected.  if averag
%		mastoids is selected this is true but if a bipolar 
%		array is selected you mayhave many more channels than the 
%		original data file.  The csdm file will however carry 
%		an updated NEeg
% 149 +NReference +NEeg +1 to  149 + NReference + NEeg +NAnalog :: any
%		analog channels
%
% Downstream analysis and visualization routines require the absolute channel
% numbers rather than the original channel numbers. 


if nargin < 2
	error('insufficient arguments');
end;

if nargin == 2
	MEG_max =  1;
	EEG_max = 1;
	EEG_reference = 'mastoids';
	add_bad_chan = [];
	max_freq = 50;
end;

if nargin == 3
	EEG_max = 1;
	EEG_reference = 'mastoids';
	add_bad_chan = [];
	max_freq = 50;
end;

if nargin == 4
	EEG_reference = 'mastoids'
	add_bad_chan = [];
	max_freq = 50;
end;

if nargin == 5
	add_bad_chan = [];
	max_freq = 50;	
end;

if nargin == 6
	max_freq = 50;
end;

if isempty(MEG_max)
	MEG_max = 1;
end

if isempty(EEG_max)
	EEG_max = 1;
end;

if isempty(EEG_reference)
	EEG_reference = 'mastoids';
end;

if isempty(add_bad_chan)
	add_bad_chan = [];
end;

if isempty(max_freq)
	max_freq = 50;
end;

if nargin < 8
	[fid,megfilename] = get_fid('rb');
	fclose(fid);
end;

meg_fid = open_file_w_byte_order(megfilename,1);

[version, NChan, NMeg, NEeg, NReference, NBad_Sensors, NTrigger, NResponse, NUtility, NAnalog, Samp_Rate, NData_Epoch, Names, NSamp, trigger, response, header_length] = rd_meg_hdr(meg_fid);

[meg_indices,meg_channels, bad_sensors,eeg_indices,eeg_channels,analog_indices,analog_channels,reference_indices,reference_channels] = parse_names(Names,NMeg,NEeg,NBad_Sensors,NAnalog,NReference);
numsamp = 4*fix(Epoch*Samp_Rate/4);
Epoch = numsamp/Samp_Rate;
NWindow = fix((NSamp-skip_samples)/numsamp);

csdmfilename = [megfilename(1:size(megfilename,2) - 4) '.csdm'];

csdm_fid = fopen(csdmfilename,'wb');
add_bad_chan = [add_bad_chan bad_sensors];
csdmversion = -1;
csdm_NMeg = 148;
%
% note: this has to be generalize dwitht he appropriate routine. 
%	to allow for bipolar pairs
csdm_NEeg = NEeg;
csdm_NReference = NReference;
csdm_NAnalog = NAnalog;

csdm_NChan = csdm_NMeg + csdm_NEeg + csdm_NReference + csdm_NAnalog;

if csdm_NChan ~= max(analog_channels);
	error('channel identification mismatch');
end;
stack = zeros(1,csdm_NChan);
csdm_NChan_file = (csdm_NChan*(csdm_NChan+1))/2;
%write out csdm header version

fwrite(csdm_fid,csdmversion,'integer*2');
avgcsdm = zeros(numsamp/4,csdm_NChan_file);
for i = 1:NWindow
	if i == 1
		trialdata = rd_meg_allch(meg_fid,header_length,NChan,numsamp,skip_samples+1);
	else
		trialdata = rd_meg_allch(meg_fid,header_length,NChan,numsamp);
	end;
	trialdata2 = trialdata(1:2:numsamp-1,:);
	trialdata_ord = fix_trial_order(trialdata2,meg_channels,meg_indices,eeg_channels,eeg_indices,analog_channels,analog_indices,reference_channels,reference_indices);
	mask = artifact_edit_meg(trialdata_ord,meg_channels,eeg_channels,MEG_max,EEG_max,add_bad_chan);
	if sum(mask) > 100
		if ~isempty(eeg_channels)
		trialdata_ord = reference_eeg_for_meg(trialdata_ord,eeg_channels,mask,EEG_reference);
		end
		avgcsdm = avgcsdm + csdm(trialdata_ord);
		stack = stack + mask;
		
	end;
end; 
icount = 1;
good_cross = zeros(1,size(avgcsdm,2));
for ichan = 1:csdm_NChan
	for jchan = ichan:csdm_NChan
		good_cross(icount) = min([stack(ichan) stack(jchan)]);
		if good_cross(icount) == 0
			good_cross(icount) = 1;
		end;
		icount = icount+1;
	end;
end;
avgcsdm = avgcsdm./(ones(size(avgcsdm,1),1)*good_cross);

NEpoch_used = max(good_cross); 

more_bad_chan = find(stack < 0.75*NEpoch_used);

add_bad_chan = [add_bad_chan more_bad_chan];

fwrite(csdm_fid,csdm_NChan_file,'integer*2');
fwrite(csdm_fid,csdm_NMeg,'integer*2');
fwrite(csdm_fid,csdm_NReference,'integer*2');
fwrite(csdm_fid,csdm_NEeg,'integer*2');
fwrite(csdm_fid,csdm_NAnalog,'integer*2');
NBad_chan = size(add_bad_chan,2);
fwrite(csdm_fid,NBad_chan,'integer*2');
fwrite(csdm_fid,add_bad_chan,'integer*2');
disp('the total set of bad channels were:')
add_bad_chan
fwrite(csdm_fid,NEpoch_used,'integer*2');
fwrite(csdm_fid,Epoch,'real*4');
max_freq = fix(max_freq*Epoch)+1;
fwrite(csdm_fid,max_freq,'integer*2');
fwrite(csdm_fid,real(avgcsdm(1:max_freq,:))','real*4');
fwrite(csdm_fid,imag(avgcsdm(1:max_freq,:))','real*4');
fclose('all');
status = 1;









