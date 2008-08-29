function [big_trial,NSamp,Samp_Rate] = rd_meg_onechan(analog_chans,megfilename);
%[analog_power,NSamp,Samp_Rate] = power_fft_analog(analog_chans,megfilename);
% rivalry script for power (analog only)
%
%analog chans = the analog channel numbers
%meg_filename = filename of MEG data, if skipped provides a gui.  
%

if nargin < 1
	error('what are the channels')
end

if nargin < 2
	[fid,megfilename] = get_fid('rb');
	fclose(fid);
end;
meg_fid = open_file_w_byte_order(megfilename,1);

[version, NChan, NMeg, NEeg, NReference, NBad_Sensors, NTrigger, NResponse, NUtility, NAnalog, Samp_Rate, NData_Epoch, Names, NSamp, trigger, response, header_length] = rd_meg_hdr(meg_fid);

[meg_indices,meg_channels, bad_sensors,eeg_indices,eeg_channels,analog_indices,analog_channels,reference_indices,reference_channels] = parse_names(Names,NMeg,NEeg,NBad_Sensors,NAnalog,NReference);

numsamp = fix(NSamp/100)
grabsamp = numsamp;
Samp_Rate = fix(Samp_Rate);
Epoch = numsamp/Samp_Rate;
NSamp = 100*numsamp;
big_trial = zeros(NSamp,size(analog_chans,2));
grab_samp = numsamp;

for itrial = 1:100
		trialdata = rd_meg_allch(meg_fid,header_length,NChan,grabsamp);
	eeg_indices = [];
	eeg_channels = [];
	trialdata_ord = fix_trial_order(trialdata,meg_channels,meg_indices,eeg_channels,eeg_indices,analog_channels,analog_indices,reference_channels,reference_indices);
	

	big_trial((itrial-1)*grabsamp+1:itrial*grabsamp,1:size(analog_chans,2)) = trialdata_ord(:,analog_chans);
end;

fclose('all')
status = 1;


















