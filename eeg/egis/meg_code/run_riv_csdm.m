function status = run_riv_csdm(skip_samples,Epoch,MEG_max,EEG_max,EEG_reference,add_bad_chan,max_freq,megfilename);
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

csdmfilename1 = [megfilename(1:size(megfilename,2) - 4) '.1.csdm'];
csdmfilename2 = [megfilename(1:size(megfilename,2) - 4) '.2.csdm'];
csdmfilenamer1 = [megfilename(1:size(megfilename,2) - 4) '.r1.csdm'];
csdmfilenamer2 = [megfilename(1:size(megfilename,2) - 4) '.r2.csdm'];
csdm_fid1 = fopen(csdmfilename1,'wb');
csdm_fid2 = fopen(csdmfilename2,'wb');
csdm_fidr1 = fopen(csdmfilenamer1,'wb');
csdm_fidr2 = fopen(csdmfilenamer2,'wb');
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
stack1 = zeros(1,csdm_NChan);
stack2 = zeros(1,csdm_NChan);
stackr1 = zeros(1,csdm_NChan);
stackr2 = zeros(1,csdm_NChan);
csdm_NChan_file = (csdm_NChan*(csdm_NChan+1))/2;
%write out csdm header version
fwrite(csdm_fid1,csdmversion,'integer*2');
fwrite(csdm_fid2,csdmversion,'integer*2');
fwrite(csdm_fidr1,csdmversion,'integer*2');
fwrite(csdm_fidr2,csdmversion,'integer*2');
avgcsdm1 = zeros(numsamp/4,csdm_NChan_file);
avgcsdm2 = zeros(numsamp/4,csdm_NChan_file);
avgcsdmr1 = zeros(numsamp/4,csdm_NChan_file);
avgcsdmr2 = zeros(numsamp/4,csdm_NChan_file);
rand_wt1 = rand(1,NWindow);
rand_wt2 = ones(1,NWindow) - rand_wt1;
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
		wt1 = size(find(response(skip_samples+(i-1)*numsamp+1:skip_samples+i*numsamp)== 1),2)/numsamp;
		wt2 = size(find(response(skip_samples+(i-1)*numsamp+1:skip_samples+i*numsamp)== 2),2)/numsamp;
		csdm_trial = csdm(trialdata_ord);
		avgcsdm1 = avgcsdm1 + wt1*csdm_trial;
		avgcsdm2 = avgcsdm2 + wt2*csdm_trial;
		avgcsdmr1 = avgcsdmr1 + rand_wt1*csdm_trial;
		avgcsdmr2 = avgcsdmr2 + rand_wt2*csdm_trial;
		stack1 = stack1 + wt1*mask;
		stack2 = stack2 + wt2*mask;
		stackr1 = stackr1 + rand_wt1*mask;
		stackr2 = stackr2 + rand_wt2*mask;
end;
end; 
icount = 1;
good_cross1 = zeros(1,size(avgcsdm1,2));
good_cross2 = zeros(1,size(avgcsdm2,2));
good_crossr1 = zeros(1,size(avgcsdm1,2));
good_crossr2 = zeros(1,size(avgcsdm1,2));
for ichan = 1:csdm_NChan
	for jchan = ichan:csdm_NChan
		good_cross1(icount) = min([stack1(ichan) stack1(jchan)]);
		if good_cross1(icount) == 0
			good_cross1(icount) = 1;
		end;
		good_cross2(icount) = min([stack2(ichan) stack2(jchan)]);
		if good_cross2(icount) == 0
			good_cross2(icount) = 1;
		end;
		good_crossr1(icount) = min([stackr1(ichan) stackr1(jchan)]);
		if good_crossr1(icount) == 0
			good_crossr1(icount) = 1;
		end;
		good_crossr2(icount) = min([stackr2(ichan) stackr2(jchan)]);
		if good_crossr2(icount) == 0
			good_crossr2(icount) = 1;
		end;
		icount = icount+1;
	end;
end;
avgcsdm1 = avgcsdm1./(ones(size(avgcsdm1,1),1)*good_cross1);
avgcsdm2 = avgcsdm2./(ones(size(avgcsdm1,1),1)*good_cross2);
avgcsdmr1 = avgcsdmr1./(ones(size(avgcsdm1,1),1)*good_crossr1);
avgcsdmr2 = avgcsdmr2./(ones(size(avgcsdm1,1),1)*good_crossr2);

NEpoch_used1 = fix(max(good_cross1)); 
NEpoch_used2 = fix(max(good_cross2)); 
NEpoch_usedr1 = fix(max(good_crossr1)); 
NEpoch_usedr2 = fix(max(good_crossr2)); 


more_bad_chan1 = find(stack1 < 0.75*NEpoch_used1);
more_bad_chan2 = find(stack2 < 0.75*NEpoch_used2);
more_bad_chanr1 = find(stackr1 < 0.75*NEpoch_usedr1);
more_bad_chanr2 = find(stackr2 < 0.75*NEpoch_usedr2);

add_bad_chan1 = [add_bad_chan more_bad_chan1];
add_bad_chan2 = [add_bad_chan more_bad_chan2];
add_bad_chanr1 = [add_bad_chan more_bad_chanr1];
add_bad_chanr2 = [add_bad_chan more_bad_chanr2];

fwrite(csdm_fid1,csdm_NChan_file,'integer*2');
fwrite(csdm_fid1,csdm_NMeg,'integer*2');
fwrite(csdm_fid1,csdm_NReference,'integer*2');
fwrite(csdm_fid1,csdm_NEeg,'integer*2');
fwrite(csdm_fid1,csdm_NAnalog,'integer*2');
NBad_chan = size(add_bad_chan1,2);
fwrite(csdm_fid1,NBad_chan,'integer*2');
fwrite(csdm_fid1,add_bad_chan1,'integer*2');
disp('the total set of bad channels were:')
add_bad_chan1
fwrite(csdm_fid1,NEpoch_used1,'integer*2');
fwrite(csdm_fid1,Epoch,'real*4');
max_freq = fix(max_freq*Epoch)+1;
fwrite(csdm_fid1,max_freq,'integer*2');
fwrite(csdm_fid1,real(avgcsdm1(1:max_freq,:))','real*4');
fwrite(csdm_fid1,imag(avgcsdm1(1:max_freq,:))','real*4');

fwrite(csdm_fid2,csdm_NChan_file,'integer*2');
fwrite(csdm_fid2,csdm_NMeg,'integer*2');
fwrite(csdm_fid2,csdm_NReference,'integer*2');
fwrite(csdm_fid2,csdm_NEeg,'integer*2');
fwrite(csdm_fid2,csdm_NAnalog,'integer*2');
NBad_chan = size(add_bad_chan2,2);
fwrite(csdm_fid2,NBad_chan,'integer*2');
fwrite(csdm_fid2,add_bad_chan2,'integer*2');
disp('the total set of bad channels were:')
add_bad_chan2
fwrite(csdm_fid2,NEpoch_used2,'integer*2');
fwrite(csdm_fid2,Epoch,'real*4');
max_freq = fix(max_freq*Epoch)+1;
fwrite(csdm_fid2,max_freq,'integer*2');
fwrite(csdm_fid2,real(avgcsdm2(1:max_freq,:))','real*4');
fwrite(csdm_fid2,imag(avgcsdm2(1:max_freq,:))','real*4');

fwrite(csdm_fidr1,csdm_NChan_file,'integer*2');
fwrite(csdm_fidr1,csdm_NMeg,'integer*2');
fwrite(csdm_fidr1,csdm_NReference,'integer*2');
fwrite(csdm_fidr1,csdm_NEeg,'integer*2');
fwrite(csdm_fidr1,csdm_NAnalog,'integer*2');
NBad_chan = size(add_bad_chanr1,2);
fwrite(csdm_fidr1,NBad_chan,'integer*2');
fwrite(csdm_fidr1,add_bad_chanr1,'integer*2');
disp('the total set of bad channels were:')
add_bad_chanr1
fwrite(csdm_fidr1,NEpoch_usedr1,'integer*2');
fwrite(csdm_fidr1,Epoch,'real*4');
max_freq = fix(max_freq*Epoch)+1;
fwrite(csdm_fidr1,max_freq,'integer*2');
fwrite(csdm_fidr1,real(avgcsdmr1(1:max_freq,:))','real*4');
fwrite(csdm_fidr1,imag(avgcsdmr1(1:max_freq,:))','real*4');

fwrite(csdm_fidr2,csdm_NChan_file,'integer*2');
fwrite(csdm_fidr2,csdm_NMeg,'integer*2');
fwrite(csdm_fidr2,csdm_NReference,'integer*2');
fwrite(csdm_fidr2,csdm_NEeg,'integer*2');
fwrite(csdm_fidr2,csdm_NAnalog,'integer*2');
NBad_chan = size(add_bad_chanr2,2);
fwrite(csdm_fidr2,NBad_chan,'integer*2');
fwrite(csdm_fidr2,add_bad_chanr2,'integer*2');
disp('the total set of bad channels were:')
add_bad_chanr2
fwrite(csdm_fidr2,NEpoch_usedr2,'integer*2');
fwrite(csdm_fidr2,Epoch,'real*4');
max_freq = fix(max_freq*Epoch)+1;
fwrite(csdm_fidr2,max_freq,'integer*2');
fwrite(csdm_fidr2,real(avgcsdmr2(1:max_freq,:))','real*4');
fwrite(csdm_fidr2,imag(avgcsdmr2(1:max_freq,:))','real*4');
fclose('all');
status = 1;










