function status = plot_egis_data(cell,obs,chans,freqs,Samp_Rate,LoPass);


if (nargin < 2);
	error('cells and observations not specified')
end;


if (nargin < 3);
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
end;
if chans == []
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
end
if nargin < 4
	freqs = [];
end

open_header_utility;
ref_trialdata = read_trial_utility(fid,Cell_DataOffset(cell),obs,NChan,Npoints(cells));
max_trialdata = max(max(abs(ref_trialdata(:,chans))));
max_trialdata = max_trialdata*1.1;
if nargin == 6
	[b,a] = butter(20,[LoPass/(Samp_Rate/2)]);
	for i = 1:size(chans,2)
		ref_trialdata = filtfilt(b,a,ref_trialdata(:,chans(i)));
	end;
end;
time = [1:size(ref_trialdata,1)];
figure
for i = 1:size(chans,2)
	hold on,plot(time,ref_trialdata(:,chans(i))+(i-1)*max_trialdata*ones(size(ref_trialdata,1),1),'w-')
	hold on, plot([min(time)-0.25 max(time)+0.25],[(i-1)*max_trialdata (i-1)*max_trialdata],'w--')
	hold on, text(max(time)+0.4,(i-1)*max_trialdata,int2str(chans(i)));
end
xlabel('Time (samples)')
ylabel('Bins')
if nargin > 3
	if freqs ~= []
		figure
		fft_trial = abs(fft(ref_trialdata(:,chans)));
		max_fft = max(max(fft_trial));
		frequency = [0:fix(size(fft_trial,1)/2)-1]/size(seconds,2);
		for i = 1:size(chans,2)
			hold on,plot(frequency,fft_trial(1:size(frequency,2),i)+(i-1)*max_fft*ones(size(frequency,2),1),'w-')
			hold on,axis([freqs(1) freqs(2) 0 size(chans,2)*max_fft ])
			hold on, plot([min(freqs)-1 max(freqs)*1.1],[(i-1)*max_fft (i-1)*max_fft],'w--')
			hold on, text(max(freqs)*0.75,(i-0.5)*max_fft,int2str(chans(i)));
			
		end
	end	
end	
	


