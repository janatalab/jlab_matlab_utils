function [names,vals] = fmri_regress_gsr(pinfo,minfo,sess)

% generates skin conducatance regressors for fmri data
% 
%   [names,vals] = fmri_regress_gsr(pinfo,minfo,sess)
% 
% takes pre-processed GSR data, resamples and normalizes it, returns it as
% a regressor
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
names = {sprintf('gsr_run%d',pinfo.irun)};
vals = [];

% init vars
TR = pinfo.TR;
nslice_per_vol = pinfo.nslice_per_vol;
actual_nvol = pinfo.actual_nvol;

% get GSR signal
gsr_epochs = pinfo.gsr_epochs;
gcol  = set_var_col_const(gsr_epochs.vars);
gsig  = gsr_epochs.data{gcol.signal}{1};
srate = gsr_epochs.data{gcol.srate}(1);

% resample, normalize GSR signal
rspl_scr = resample(double(gsig),nslice_per_vol/TR,srate);
if length(rspl_scr) > actual_nvol
  rspl_scr(actual_nvol+1:end) = [];
elseif length(rspl_scr) < actual_nvol
  rspl_scr(end+1:actual_nvol) = 0;
end
vals = zscore(rspl_scr);
