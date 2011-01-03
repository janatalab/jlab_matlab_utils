function [names,vals] = fmri_regress_hrflux(pinfo,minfo,sess)

% generates heart rate flucuation regressors for fmri data
% 
%   [names,vals] = fmri_regress_hrflux(pinfo,minfo,sess)
% 
% per Shmueli et al (2007) NeuroImage
% 
% heart rate regressors (representing low-frequency fluctuation in
% heart rate during an imaging run) were found to explain significant
% variance within a volume-based GLM. we include an adaptation of the
% Shmueli method here.
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array of six regressor names
%   vals = volume X regressor matrix containing motion regressors
% 
% FB 2009.11.05

% init output vars
names = {sprintf('hrv_run%d',pinfo.irun)};
vals = zeros(pinfo.scanner.actual_nvol,1);

% init vars
TR = pinfo.scanner.TR;
nvol = pinfo.scanner.actual_nvol;

% get cardiac data
card_epochs = pinfo.card_epochs;
ccol = set_var_col_const(card_epochs.vars);
pidxs = card_epochs.data{ccol.peakidxs}{1};
srate = card_epochs.data{ccol.srate}(1);

peaks = pidxs/srate; % peak times in seconds
hr_vect = [0; 1./diff(peaks)]; % peak-to-peak heart rate

% volume onset times, in seconds
vol_onsets = 0:TR:(nvol-1)*TR;

% per Shmueli 2007, choose the closest momentary HR to the beginning
% of each volume
for ii=1:nvol
  c1 = find(peaks <= vol_onsets(ii),1,'last');
  c2 = find(peaks >= vol_onsets(ii),1,'first');
  if isempty(c1) || c1 == 0
    c1 = 1;
  elseif c1==length(peaks)
    c1 = length(peaks)-1;
  end
  if isempty(c2)
    c2 = c1 + 1;
  elseif c2==1
    c2 = 2;
  elseif c2 > length(peaks)
    c2 = length(peaks);
  end
  if abs(peaks(c1)-vol_onsets(ii)) < abs(peaks(c2)-vol_onsets(ii))
    vals(ii) = hr_vect(c1);
  else
    vals(ii) = hr_vect(c2);
  end
end
