function avg_coherence = eeg_coherence(ravgcsdm,iavgcsdm,bad_chan);
%avg_coherence = eeg_coherence(ravgcsdm,iavgcsdm,bad_chan);
%
%avgcsdm= average cross spectral density matrix
%bad_chan = bad_channels (optional)
%avg_coherence = coherence matrix
%
if nargin < 2
	error('you missing either the real or the imaginary part');
end;
if nargin == 2
	bad_chan = [];
end;
if size(ravgcsdm,2) == 8385
	NChan = 129;
else
	error('unknown number of channels');
end;
mask = ones(1,NChan);
mask(bad_chan) = zeros(1,size(bad_chan,2));
good_chan = find(mask);
ch_pair_indices;

avg_coherence = zeros(size(ravgcsdm,1),size(ravgcsdm,2));

for i = 1:size(good_chan,2)
	for j = i:size(good_chan,2)
		avg_coherence(:,chpair(good_chan(i),good_chan(j)))= ((ravgcsdm(:,chpair(good_chan(i),good_chan(j))).^2)+(iavgcsdm(:,chpair(good_chan(i),good_chan(j))).^2))./ravgcsdm(:,chpair(good_chan(i),good_chan(i)))./ravgcsdm(:,chpair(good_chan(j),good_chan(j)));
	end;
end;

avg_coherence(find(isnan(avg_coherence))) = zeros(size(find(isnan(avg_coherence))));

