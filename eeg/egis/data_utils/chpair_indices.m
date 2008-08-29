function chpair = chpair_indices(nchan)
% chpair = chpair_indices(nchan);
%
%

if nargin < 1
  nchan = 129;
end

icount = 1;
for ichan = 1:nchan
  for jchan = ichan:nchan
    chpair(ichan, jchan) = icount;
    chpair(jchan,ichan) = icount;
    icount = icount + 1;
  end
end
