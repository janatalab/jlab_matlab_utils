function [cross_trial, fft_trial] = fastcsdm(trialdata)
% [cross_trial, fft_trial] = fastcsdm(trialdata);
%
% This is faster than csdm up to 256 pt FFTs.  After that it actually seems to
% get slower.  4 fold speed increase on 64pt FFTs
%
  
trialdata = zeromean(trialdata);
fft_trial = fft(trialdata)./(size(trialdata,1)+1);

nfreq = size(trialdata,1)/2;
nchan = size(trialdata,2);

cross_trial = zeros(nfreq,(nchan.^2+nchan)/2); 
cross_prod = zeros(nchan);

relevant_indices = find(tril((1:nchan)'*(1:nchan)));

if strcmp(computer, 'MAC2')
  disp('Run this on a real machine')
else
  for ifreq = 1:nfreq
    chanvect = fft_trial(ifreq,:);
    cross_prod = chanvect.'*conj(chanvect);
    cross_trial(ifreq,:) = cross_prod(relevant_indices)';
  end;
end
