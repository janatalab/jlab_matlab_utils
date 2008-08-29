function pairlist =  get_all_pairs2(chans,nchan);

chpair = channel_pairs(nchan);
ic = 1;
for i = 1:size(chans,2)-1
	for j = i+1:size(chans,2)
		pairlist(ic) = chpair(chans(i),chans(j));
		ic = ic+1;
	end;
end;

	
