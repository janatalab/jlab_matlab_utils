function coherence_array = get_eeg_coherence(fid,freqs,obs,nfiles,NChan,NFreq,Epoch)
%
%
if fid < 0 
	error('what the hell are you doing')
end;
obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
icount = 1;
for i = 1:nfiles
	if obs_mask(i) > 0
	avgcoherence = fread(fid,[NChan(i),NFreq(i)],'float');
	avgcoherence = avgcoherence';
	freq = freqs*Epoch(i) + ones(1,size(freqs,2));
	coherence_array(icount:icount+size(freqs,2)-1,:) = avgcoherence(freq,:);
	icount = icount+size(freqs,2);
	else
	skip_bytes = NChan(i)*NFreq(i)*4;
	fseek(fid,skip_bytes,'cof');
	end
end;
if (icount-1 ~= size(freqs,2)*size(obs,2))
	error('coherence read had a problem')
end; 
frewind(fid);
