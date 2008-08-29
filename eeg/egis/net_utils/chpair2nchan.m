function nchan = chpair2nchan(nchpair);

ichan = [1:256];
chchan = (ichan.*(ichan+ones(1,256)))/2;
nchan = ichan(find(chchan == nchpair));