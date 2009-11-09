function rvals = fmri_convolve_regress(ons, dur, amp, TR, dt, k)

% convolves regressors with SPM canonical HRF
%
%   rvals = frmi_convolve_regress(ons,dur,amp,TR,dt,k)
% 
% The core of this code is copied from
% spm_get_ons.m in the SPM5 distribution
% Assumes onsets and durations are given in seconds  
% dt=1/16 s by default
% 
% 2009.11.05 FB - adapted from fmri_spm_generate_regress, w/c was adapted
% from code written by PJ for autobio, which was adapted from spm_get_ons.m

TR = 1/TR;  % TR - this needs to be pulled from the protocol info
dt = 1/dt;  % This needs to be pulled from the fMRI specs
T = 1/dt;
  
% create stimulus functions (32 bin offset)
%===============================================================
ton       = round(ons*TR/dt) + 32;			% onsets
tof       = round(dur*TR/dt) + ton + 1;			% offset
sf        = sparse((k*T + 128),1);
ton       = max(ton,1);
tof       = max(tof,1);
for j = 1:length(ton)
  if numel(sf)>ton(j),
    sf(ton(j),:) = sf(ton(j),:) + amp(j);
  end;
  if numel(sf)>tof(j),
    sf(tof(j),:) = sf(tof(j),:) - amp(j);
  end;
end
sf        = cumsum(sf);					% integrate
sf        = sf(1:(k*T + 32),:);				% stimulus

% Convolve with the desired basis function
rvals = conv(full(sf), spm_hrf(dt/TR));

% Resample onto timescale of scans
rvals = rvals([0:(k-1)]*T+1+32);
