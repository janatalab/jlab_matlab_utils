function [names,vals] = fmri_regress_response(pinfo,minfo,sess)

% generates response onset regressors for fmri data
% 
%   [names,vals] = fmri_regress_response(pinfo,minfo,sess)
% 
% regressor to model button pressing
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
names = {'response'};
vals = [];

pc = set_var_col_const(pinfo.vars);
respidxs = find(~cellfun(@isempty,pinfo.data{pc.RESP_TIME}));
resp_onsets = [];
for ievt = 1:length(respidxs)
  evt = pinfo.data{pc.RESP_TIME}(respidxs(ievt));
  evt_onset = pinfo.data{pc.RUN_REL_TIME}(respidxs(ievt));
  for ion=1:length(evt{1})
    onset = evt{1}(ion) + evt_onset;
    resp_onsets = [resp_onsets; onset];
  end
end
onsets = resp_onsets/1000;

% Durations and amplitude
durs = zeros(size(onsets));
amp  = ones(size(onsets));

% Now build the regressor
vals = fmri_convolve_regress(onsets,durs,amp,pinfo.TR,pinfo.dt,pinfo.actual_nvol);
