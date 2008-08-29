function status = plot_meg_data(samples, chans, freqs, Lo_Pass,Hi_Pass,megfname);

%plot_meg_data(samples, chans, freqs,Lo_Pass,Hi_Pass,
%              megfname);
%plots meg data and its spectra for a given range of samples and channels. 
%artifact editing on the basis of amplitude criterion in pico (10^-12) Tesla
%can be used to blank out bad channels or tune editing parameters.
%
%samples = [min_samp max_samp]
%chans = [chan1 chan2 chan3]
%freqs = [freq_min freq_max]
%Lo_Pass = low_pass filter frequency (defaults to 50)
%Hi_Pass = hi_pass filter frequency (defaults to 0)
%megfname = meg_filename (if not provided a gui comes up)
%
if nargin < 1
	error('samples not specified');
end;

if nargin == 1
	chans = 1;
	freqs = [0 50];
	Lo_Pass = 50;
	Hi_Pass = 1;
end;
if nargin == 2
	freqs = [0 50];
	Lo_Pass = 50;
	Hi_Pass = 1;
end;
if nargin == 3
	Lo_Pass = 50;
	Hi_Pass = 1;
end;
if nargin == 4
	Hi_Pass = 1;
end;

if nargin < 6
	[meg_fid, megfname]  = get_fid('rb');
	fclose(meg_fid);
end;

meg_fid = open_file_w_byte_order(megfname,1);

if isempty(Hi_Pass);
	Hi_Pass= 1;
end;
if isempty(Lo_Pass)
	Lo_Pass = 50;
end;
if isempty(freqs);
	freqs = [0 40];
end

if isempty(chans)
	chans = 1;
end;
[version, NChan, NMeg, NEeg, NReference, NBad_Sensors, NTrigger, NResponse, NUtility, NAnalog, Samp_Rate, NData_Epoch, Names, NSamp, trigger, response, header_length] = rd_meg_hdr(meg_fid);
[meg_indices,meg_channels, bad_sensors,eeg_indices,eeg_channels,analog_indices,analog_channels,reference_indices,reference_channels] = parse_names(Names,NMeg,NEeg,NBad_Sensors,NAnalog,NReference);
trialdata = rd_meg_allch(meg_fid,header_length,NChan,samples(2) - samples(1)+1,samples(1));
trialdata_ord = fix_trial_order(trialdata,meg_channels,meg_indices,eeg_channels,eeg_indices,analog_channels,analog_indices,reference_channels,reference_indices);
trialdata_ord = zeromean(trialdata_ord);
max_trialdata = max(max(abs(trialdata_ord(:,chans))));
max_trialdata = max_trialdata*1.1;
step = 1/Samp_Rate;
time = [step*samples(1):step:step*samples(2)];
[blow,alow] = butter(10,Lo_Pass/(fix(Samp_Rate)/2));
[bhigh,ahigh] = butter(2,[Hi_Pass]/(fix(Samp_Rate)/2),'high');
for i = 1:size(chans,2)
	trialdata_ord(:,chans(i)) = filtfilt(blow,alow,trialdata_ord(:,chans(i)));
	trialdata_ord(:,chans(i)) = filtfilt(blow,alow,trialdata_ord(:,chans(i)));
end;
figure
for i = 1:size(chans,2)
	hold on,plot(time,trialdata_ord(:,chans(i))+(i-1)*max_trialdata*ones(size(trialdata_ord,1),1),'k-')
	hold on, plot([min(time)-0.25 max(time)+0.25],[(i-1)*max_trialdata (i-1)*max_trialdata],'k--')
	hold on, text(max(time)+0.4,(i-1)*max_trialdata,int2str(chans(i)));
end
total_samp = samples(2)-samples(1)+1;
xlabel('Time (seconds)')
ylabel('PicoTesla')
	if ~isempty(freqs)
		figure
		if freqs(1) == 0
			freqs(1) = 0.1;
		end;
		fft_trial = abs(fft(trialdata_ord(:,chans)))/(total_samp+1);
		binmin = ceil(freqs(1)*(total_samp/(fix(Samp_Rate))));
		binmax = fix(freqs(2)*(total_samp/(fix(Samp_Rate))));

		max_fft = max(max(fft_trial(binmin:binmax,:)));
		frequency = [0:fix(size(fft_trial,1)/2)-1]/(total_samp/fix(Samp_Rate));
		for i = 1:size(chans,2)
			hold on,plot(frequency,fft_trial(1:size(frequency,2),i)+(i-1)*max_fft*ones(size(frequency,2),1),'k-')
			hold on,axis([freqs(1) freqs(2) 0 size(chans,2)*max_fft ])
			hold on, plot([min(freqs)-1 max(freqs)*1.1],[(i-1)*max_fft (i-1)*max_fft],'k--')
			hold on, text(max(freqs)*0.75,(i-0.5)*max_fft,int2str(chans(i)));
						
		end
		xlabel('Frequency(Hz)')
		ylabel('Amplitude');
	end	
	
	





