% ch_pair_indices.m

icount = 1;
for ichan = 1:129
  for jchan = ichan:129
    chpair(ichan, jchan) = icount;
	chpair(jchan,ichan) = icount;
    icount = icount + 1;
  end
end
