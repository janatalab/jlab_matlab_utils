function status = plot_eeg_coherence(freqs,obs,chans,lines,coherencefname);
%status = plot_eeg_coherence(freqs,obs,chans,lines,coherencefname);
%
%freqs = frequencies to plot, e.g. [3 30] limits the plot to 3 to 30 Hz  
%obs = observations in the power file to plot (defaults to all)
%chans = channel pairs to plot (defaults to all pairs of F3 F4 P3 P4);
%	e.g. [26 124; 32 34;], i.e. the pairs are separated by a semicolon
% nobs <= 5. 
%lines = 'black' or 'color'
%coherencefname = coherencefilename (can be gui by skipping or passing a blank)
%
if nargin < 5
	[fid] = get_fid('rb');
else
	fid = fopen(coherencefname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -2
	error('this is not a coherence file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [0 40];
	obs = [1:nfiles];
	chans = [25 124;25 62;25 86;124 62;124 86;62 86];
	lines = 'black';
end;
if nargin == 1
	obs = [1:nfiles];
	chans = [25 124;25 62;25 86;124 62;124 86;62 86];
	lines = 'black';
end;		
if nargin == 2
	chans = [25 124;25 62;25 86;124 62;124 86;62 86];
	lines = 'black';
end;	
if nargin == 3
	lines = 'black';
end;

if isempty(lines)
	lines = 'black';
end
if isempty(chans)
	chans = [25 124;25 62;25 86;124 62;124 86;62 86];
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
if size(chans,1) > 1
layout = ceil(sqrt(size(chans,1)+1));
else
layout = 1;
end;
ch_pair_indices;
if strcmp(lines,'black')
line_type = ['k- ';'k--';'k-.';'kx ';'ko '];
else	
line_type = ['k-';'r-';'b-';'g-';'m-';'c-'];
end;
figure
coherencemax = 1;
obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
iline = 1;
for i = 1:nfiles
	if obs_mask(i) > 0
	avgcoherence = fread(fid,[NChan(i),NFreq(i)],'float');
	avgcoherence = avgcoherence';
	binmin = freqs(1)*Epoch(i) + 1;
	binmax = freqs(2)*Epoch(i) + 1;
	freqstep = 1/Epoch(i);
	freq = [freqs(1):freqstep:freqs(2)];
	for ic = 1:size(chans,1) 
		hold on, subplot(layout,layout,ic), plot(freq,avgcoherence(binmin:binmax,chpair(chans(ic,1),chans(ic,2))),line_type(iline,:))
		if layout < 3
			title(['Channels ' int2str(chans(ic,1)) ':' int2str(chans(ic,2))])
			xlabel('Frequency(Hz)')
			ylabel('Coherence')
		end	
	end
	iline = iline + 1;
	else
	skip_bytes = NChan(i)*NFreq(i)*4;
	fseek(fid,skip_bytes,'cof');
	end;
end	
for ic = 1:size(chans,1)
	subplot(layout,layout,ic), axis([0 freqs(2) 0 1]);
end;
maxcoherence = 1;
if layout < 2
	for i = 1:size(obs,2)
		hold on, plot([0.6*freqs(2) 0.7*freqs(2)],[(1-0.1*i)*maxcoherence  (1-0.1*i)*maxcoherence],line_type(i,:))
		text(0.75*freqs(2),(1-0.1*i)*maxcoherence,deblank(setstr(obs_labels(obs(i),:))));
	end;
end 
if layout > 2
	for i = 1:size(chans,1)
		hold on, subplot(layout,layout,i), text(0.6*freqs(2),0.6*maxcoherence,[int2str(chans(i,1)) ' : ' int2str(chans(i,2))]);
	end;
status = 1;
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

