% ch_pair_indices.m

icount = 1;
for ichan = 1:65
  for jchan = ichan:65
    chpair(ichan, jchan) = icount;
	chpair(jchan,ichan) = icount;
    icount = icount + 1;
  end
end
