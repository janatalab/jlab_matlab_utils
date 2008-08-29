function avg_coherence = meg_coherence(ravgcsdm,iavgcsdm,bad_chan);
NChan = 148;
	mask = ones(1,NChan);
	mask(bad_chan) = zeros(1,size(bad_chan,2));
	good_chan = find(mask);
	chpair = channel_pairs(NChan);
	avg_coherence = zeros(size(ravgcsdm,1),size(ravgcsdm,2));

	for i = 1:size(good_chan,2)
	for j = i:size(good_chan,2)
		avg_coherence(:,chpair(good_chan(i),good_chan(j)))= ((ravgcsdm(:,chpair(good_chan(i),good_chan(j))).^2)+(iavgcsdm(:,chpair(good_chan(i),good_chan(j))).^2))./ravgcsdm(:,chpair(good_chan(i),good_chan(i)))./ravgcsdm(:,chpair(good_chan(j),good_chan(j)));
	end;
end;

	avg_coherence(find(isnan(avg_coherence))) = zeros(size(find(isnan(avg_coherence))));
