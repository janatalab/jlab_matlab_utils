function cohmat = coh_from_csdm(csdmdata)
%  cohmat = coh_from_csdm(csdmdata)
%
%  csdmdata, dim 1 = freqbands
%  csdmdata, dim 2 = crossed chans
%  csdmdata, dim 3 (optional) = cells

% 01/30/00 PJ

nfreq = size(csdmdata, 1);

if (size(csdmdata,2) ~= (129^2 + 129)/2)
  error('Not equiped to handle non 129 channel data')
end

ch_pair_indices
power_indices = diag(chpair);

lower_indices = find(tril(ones(129)));

cohmat = zeros(size(csdmdata));

ncells = size(csdmdata,3);

for icell = 1:ncells
  for ifreq = 1:nfreq
    datavect = csdmdata(ifreq,:,icell);
    powervect = csdmdata(ifreq,power_indices,icell);
    
    cross_power = powervect' * powervect;
    
    cohmat(ifreq,:,icell) = (datavect .* conj(datavect))./cross_power(lower_indices)';
  end
end