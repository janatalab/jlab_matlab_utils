function status = plot_meg_coherence(freqs,chans,obs,lines,coherencefilename);
%
% plots meg coherence
%
% freqs = [min_freq mex_freq];
% chans = channel pair list (absolute)
% obs = observation number
% lines = 'color' or 'black'
% coherencefilename = can be gui
%

if nargin < 3
	error('insufficient input arguments');
end;

if nargin < 5
	[fid,coherencefilename] = get_fid('rb');
	fclose(fid);
end;

if nargin == 3
	lines = 'black';
end;

if isempty(lines)
	lines = 'black';
end

if size(obs,2) > 4
	error('too many observations for a single plot');
end;
if lines  == 'black'
	line_type = ['k- ';'k--';'k-.';'k- '];
else
	line_type = ['k-';'r-';'b-';'g-'];
end;
pfid = open_file_w_byte_order(coherencefilename,-3);
if size(chans,1) > 1
layout = ceil(sqrt(size(chans,1)+1));
else
layout =1;
end
[pversion,nfile,NChan,NMeg,NReference,NEeg,NAnalog,NBad_chan,bad_chan,NEpoch,Epoch,nfreq,obs_labels] = rd_meg_anal_hdr(pfid);
obs_mask = zeros(1,nfile);
all_NChan = NChan.*(NChan+ones(1,size(NChan,2)))/2;
obs_mask(obs) = ones(1,size(obs,2));
coherencemax = 1;
figure
iline = 1;
for i = 1:nfile
	if obs_mask(i) > 0
	avgcoherence = fread(fid,[all_NChan(i),nfreq(i)],'float');
	avgcoherence = avgcoherence';
	binmin = fix(freqs(1)*Epoch(i) + 1);
	binmax = fix(freqs(2)*Epoch(i) + 1);
	freqstep = 1/Epoch(i);
	freq = [(binmin-1):(binmax-1)]*(1/(Epoch(i)));
	chpair = channel_pairs(NChan(i));
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
	skip_bytes = NChan(i)*nfreq(i)*4;
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



