function status = plot_eeg_power(freqs,obs,chans,lines,power_type,powerfname);
%status = plot_eeg_power(freqs,obs,chans,lines,power_type,powerfname);
%
%freqs = frequencies to plot, e.g. [3 30] limits the plot to 3 to 30 Hz  
%obs = observations in the power file to plot (defaults to all)
%chans = channels to plot (defaults to 10_20);
% nobs <= 5. 
%lines = 'black' or 'color'
%power_type = can be either 'power', 'amplitude', or 'logpower'
%powerfname = powerfilename (can be gui by skipping or passing a blank)
%
if nargin < 6
	[fid] = get_fid('rb');
else
	fid = fopen(powerfname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -3
	error('this is not a power file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [0 40];
	obs = [1:nfiles];
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
	power_type = 'amplitude';
	lines = 'black';
end;
if nargin == 1
	obs = [1:nfiles];
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
	power_type = 'amplitude';
	lines = 'black';
end;		
if nargin == 2
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
	power_type = 'amplitude';
	lines = 'black';
end;	
if nargin == 3
	power_type = 'amplitude';
	lines = 'black';
end;
if nargin == 4
	power_type = 'amplitude';
end;
if isempty(lines)
	lines = 'black';
end
if isempty(power_type)
	power_type = 'amplitude';
end;
if isempty(chans)
	
	chans = [9 11 23 25 34 37 46 58 60 62 71 84 86 97 105 109 122 124 129];
end
if isempty(obs)
	obs = [1:nfiles];
end;
if isempty(freqs)
	freqs = [0 40];
end;
if size(obs,2) > 5
	fclose('all')
	error('too many observations > 5 to be plotted in a single plot')
end;
if size(chans,2) > 1
layout = ceil(sqrt(size(chans,2)+1));
else
layout =1;
end
if strcmp(lines,'black')
line_type = ['k- ';'k--';'k-.';'kx ';'ko '];
else	
line_type = ['k-';'r-';'b-';'g-';'m-';'c-'];
end;
figure
powermax = 0;
obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
maxpower = 0;
iline = 1;
for i = 1:nfiles
	power = fread(fid,[NChan(i),NFreq(i)],'float');
	if obs_mask(i) > 0
	power = power';
	if strcmp(power_type,'amplitude')
		power = sqrt(power);
	end;
	binmin = freqs(1)*Epoch(i) + 1;
	binmax = freqs(2)*Epoch(i) + 1;
	freqstep = 1/Epoch(i);
	freq = [freqs(1):freqstep:freqs(2)];
	for ic = 1:size(chans,2) 
		if strcmp(power_type,'logpower')
			hold on,subplot(layout,layout,ic),  semilogy(freq,power(binmin:binmax,chans(ic)),line_type(iline,:))
			if layout < 3
			title(['Channel ' int2str(chans(ic))])
			xlabel('Frequency(Hz)')
			ylabel('Power')
			end;
		else
			hold on, subplot(layout,layout,ic), plot(freq,power(binmin:binmax,chans(ic)),line_type(iline,:))
			if layout < 3
				title(['Channel ' int2str(chans(ic))])
				xlabel('Frequency(Hz)')
				if strcmp(power_type,'amplitude')
					ylabel('Amplitude');
				else
					ylabel('Power')
				end				
			end;
		end	
	end
	mpower = max(max(power(binmin:binmax,chans)));
	maxpower = max([maxpower mpower]);
	iline = iline + 1;
	end;
end	
for ic = 1:size(chans,2)
	subplot(layout,layout,ic), axis([0 freqs(2) 0 1.05*maxpower]);
end;
if layout < 2
	for i = 1:size(obs,2)
		hold on, plot([0.6*freqs(2) 0.7*freqs(2)],[(1-0.1*i)*maxpower  (1-0.1*i)*maxpower],line_type(i,:))
		text(0.75*freqs(2),(1-0.1*i)*maxpower,deblank(setstr(obs_labels(obs(i),:))));
	end;
end
if layout > 2
	for i = 1:size(chans,2)
		hold on, subplot(layout,layout,i), text(0.8*freqs(2),0.8*maxpower,int2str(chans(i)));
	end;
end;
if layout > 1
	for i = 1:size(obs,2)
		hold on, subplot(layout,layout,layout*layout),plot([0.2 1],[1 1]*(1.1-0.2*i),line_type(i,:)), axis([0 1.5 0 1])
		hold on, subplot(layout,layout,layout*layout),text(1.2,(1.1-0.2*i),deblank(setstr(obs_labels(obs(i),:))))
		hold on, subplot(layout,layout,layout*layout),axis('off')
	end;
	hold on, subplot(layout,layout,layout*layout),axis('off');
end;
status = 1;


		



