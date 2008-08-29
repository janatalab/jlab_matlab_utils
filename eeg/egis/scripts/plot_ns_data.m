function bad_chan= plot_ns_data(seconds,reference,chans,Max_microv,freqs,rawfname);
%bad_chan= plot_ns_data(seconds,reference,channels,Max_microv,freqs,rawfname);
%
%seconds = seconds from start of file to grab data, e.g., [1:4] or [5] or [97:111]
%reference = reference to use (see new_reference.m) (optional)
%		options: 'average', 'avgmast', 'perimeter', 'vertex',or a channel list
% 		defaults to 'average'
%		can also be used to invoke laplacian: 'laplacian'
%		can also be used to invoke cortical filtering: 'cortical05'
%			or 'cortical10' where the number indicates the
%			smoothing level
%data is alway zeromeaned channel by channel
%channels = list of channels to plot (defaults to 10_20)
%Max_microv = maximum microvolt level at channels (optional)
%				defaults to 100
%freqs = frequency plot if skipped or [] then 
%no spectral plot
%rawfname = filename (optional)
%
if nargin < 1
	error('seconds not specified');
end;
if (nargin < 2);
	reference = 'average';
end;
if isempty(reference);
	reference = 'average';
end;
if (nargin < 3);
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
end;
if chans == []
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
end
if (nargin < 4)
	Max_microv = 100;
end
if Max_microv == [];
	Max_microv = 100;
end;
if nargin < 6
	[ref_trialdata, bad_chan] = grab_ns_data(seconds,reference,Max_microv);
else
	[ref_trialdata, bad_chan] = grab_ns_data(seconds,reference,Max_microv,rawfname);
end;
max_trialdata = max(max(abs(ref_trialdata(:,chans))));
max_trialdata = max_trialdata*1.1;
step = size(seconds,2)/size(ref_trialdata,1);
time = [step:step:size(seconds,2)];
time = time+min(seconds)*ones(1,size(time,2));
figure
for i = 1:size(chans,2)
	hold on,plot(time,ref_trialdata(:,chans(i))+(i-1)*max_trialdata*ones(size(ref_trialdata,1),1),'k-')
	hold on, plot([min(time)-0.25 max(time)+0.25],[(i-1)*max_trialdata (i-1)*max_trialdata],'k--')
	hold on, text(max(time)+0.4,(i-1)*max_trialdata,int2str(chans(i)));
end
xlabel('Time (seconds)')
ylabel('Microvolts')
if nargin > 4
	if freqs ~= []
		figure
		fft_trial = abs(fft(ref_trialdata(:,chans)));
		max_fft = max(max(fft_trial));
		frequency = [0:fix(size(fft_trial,1)/2)-1]/size(seconds,2);
		for i = 1:size(chans,2)
			hold on,plot(frequency,fft_trial(1:size(frequency,2),i)+(i-1)*max_fft*ones(size(frequency,2),1),'k-')
			hold on,axis([freqs(1) freqs(2) 0 size(chans,2)*max_fft ])
			hold on, plot([min(freqs)-1 max(freqs)*1.1],[(i-1)*max_fft (i-1)*max_fft],'k--')
			hold on, text(max(freqs)*0.75,(i-0.5)*max_fft,int2str(chans(i)));
			
		end
	end	
end	
	

