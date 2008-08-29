function chpair = channel_pairs(NChan);


icount = 1;
for ichan = 1:NChan
  for jchan = ichan:NChan
    chpair(ichan, jchan) = icount;
	chpair(jchan,ichan) = icount;
    icount = icount + 1;
  end
end
