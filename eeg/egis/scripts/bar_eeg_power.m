function status = bar_eeg_power(freqs,obs,chans,base_obs,lines,power_type,powerfname);
%status =bar_eeg_power(freqs,obs,chans,base_obs,lines,power_type,powerfname)
%
%freqs = list of freqs, e.g. [3 30] limits the plot to 3 to 30 Hz  
%obs = observations in the power file corresponding to freqs
%note: nfreqs == nobs
%chans = channels to plot; channel list separatedby a semicolon will be
%averaged together .  standard array names can be used.  
%base_obs = observation to plot at each freq under freqs asa line graph.  
%lines = 'black' or 'color'
%power_type = can be either 'power', 'amplitude', or 'logpower'
%powerfname = powerfilename (can be gui by skipping or passing a blank)
%

if nargin < 6
	error('insufficient number of arguments')
end
if size(obs,2) ~= size(freqs,2) 
	error('frequencies and observations dont match')
end
if nargin < 7
	[fid] = get_fid('rb');

else
	fid = fopen(powerfname,'rb');
end;
version = fread(fid,1,'int16');
if version ~= -3
	error('this is not a power file');
end;

 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);

layout = ceil(sqrt(size(chans,2)));

if strcmp(lines,'black')
line_type1 = ['w-'];
line_type2 = ['k--'];
else	
line_type1 = ['r-'];
line_type2 = ['b--'];
end;

figure

obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
maxpower = 0;
plot_power = zeros(size(obs,2),size(chans,2));
iline = 1;
icount = 1;
for i = 1:nfiles
	power = fread(fid,[NChan(i),NFreq(i)],'float');
	if strcmp(power_type,'amplitude')
			power = sqrt(power);
	end;
	if i == base_obs
		power_base = power';
	end;
	if obs_mask(i) > 0
		power = power';
		
		plot_power(icount,1:size(chans,2)) = power(freqs(icount)*Epoch(i)+1,chans);
		icount = icount + 1;
	end;
max_power = max(max(plot_power));
end;
icount = 1;
for i = 1:size(chans,2)
	subplot(layout,layout,icount), bar(freqs,plot_power(:,icount),line_type1)
	hold on, plot(freqs,power_base(freqs*Epoch(base_obs)+ones(1,size(freqs,2)),chans(icount)),line_type2);
	hold on, axis([min(freqs)-1 max(freqs)+1 0 max_power]);
	title(['Channel ' int2str(chans(icount))]);
	icount = icount + 1;
	
end;
fclose('all')
status = 1;
	


		







