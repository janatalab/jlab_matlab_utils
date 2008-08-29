function status = plot_meg_power(freqs,chans,obs,power_type,lines,powerfilename);
%
% plots meg power
%
% freqs = [min_freq mex_freq];
% chans = channel list (absolute)
% obs = observation number
% power_type = 'amplitude','power','logpower'
% lines = 'color' or 'black'
% powerfilename = can be gui
%

if nargin < 4
	error('insufficient input arguments');
end;

if nargin < 6
	[fid,powerfilename] = get_fid('rb');
	fclose(fid);
end;

if nargin == 4
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
pfid = open_file_w_byte_order(powerfilename,-2);
if size(chans,2) > 1
lay_out = ceil(sqrt(size(chans,2)+1));
else
lay_out =1;
end
[pversion,nfile,NChan,NMeg,NReference,NEeg,NAnalog,NBad_chan,bad_chan,NEpoch,Epoch,nfreq,obs_labels] = rd_meg_anal_hdr(pfid);
obs_mask = zeros(1,nfile);
obs_mask(obs) = ones(1,size(obs,2));
maxpower = 0;
figure
iline = 1;
for i = 1:nfile
	power = fread(pfid,[NChan(i),nfreq(i)],'real*4');
	if obs_mask(i) > 0
	power = power';
	if strcmp(power_type,'amplitude')
		power = sqrt(power);
	end;
	binmin = fix(freqs(1)*Epoch(i) + 1);
	binmax = fix(freqs(2)*Epoch(i) + 1);
	freqstep = 1/Epoch(i);
	freq = [(binmin-1):(binmax-1)]*(1/(Epoch(i)));
	for ic = 1:size(chans,2) 
		if strcmp(power_type,'logpower')
			hold on,subplot(lay_out,lay_out,ic),  semilogy(freq,power(binmin:binmax,chans(ic)),line_type(iline,:))
			if lay_out < 3
			if chans(ic) <= 148
			title(['MEG ' int2str(chans(ic))])
			elseif chans(ic) <= 148+NReference
			title(['Ref ' int2str(chans(ic))])
			elseif chans(ic) <= 148+NReference+NEeg
			title(['EEG ' int2str(chans(ic))])
			else
			title(['Ana ' int2str(chans(ic))])
			end;
			xlabel('Frequency(Hz)')
			ylabel('Power')
			end;
		else
			hold on, subplot(lay_out,lay_out,ic), plot(freq,power(binmin:binmax,chans(ic)),line_type(iline,:))
			if lay_out < 3
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
	subplot(lay_out,lay_out,ic), axis([freqs(1) freqs(2) 0 1.05*maxpower]);
end;
if lay_out < 2
	for i = 1:size(obs,2)
		hold on, plot([0.6*freqs(2) 0.7*freqs(2)],[(1-0.1*i)*maxpower  (1-0.1*i)*maxpower],line_type(i,:))
		
		text(0.75*freqs(2),(1-0.1*i)*maxpower,obs_labels(obs(i),:));
	end;
end
if lay_out > 2
	for i = 1:size(chans,2)
		if chans(i) <= 148
		chan_label = ['MEG ' int2str(chans(i))];
		elseif chans(i) <= 148 + NReference
		chan_label = ['Ref ' int2str(chans(i))];
		elseif chans(i) <= 148 + NReference +NEeg
		chan_label = ['EEG ' int2str(chans(i))];
		else
		chan_label = ['Ana ' int2str(chans(i))];
		end 
		hold on, subplot(lay_out,lay_out,i), text(0.8*freqs(2),0.8*maxpower,chan_label);
	end;
end;
if lay_out > 1
	for i = 1:size(obs,2)
		hold on, subplot(lay_out,lay_out,lay_out*lay_out),plot([0.2 1],[1 1]*(1.1-0.2*i),line_type(i,:)), axis([0 1.5 0 1])
		hold on, subplot(lay_out,lay_out,lay_out*lay_out),text(1.2,(1.1-0.2*i),obs_labels(obs(i),:))
		hold on, subplot(lay_out,lay_out,lay_out*lay_out),axis('off')
	end;
	hold on, subplot(lay_out,lay_out,lay_out*lay_out),axis('off');
end;
status = 1;





