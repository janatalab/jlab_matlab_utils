function [ons,durs] = fmri_cond_response(pinfo,minfo,sess)

% generates response onset conditions for fmri data
% 
%   [ons,durs] = fmri_cond_response(pinfo,minfo,sess)
% 
% condition generator to model button pressing
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
ons = resp_onsets/1000;

% Durations and amplitude
durs = zeros(size(ons));
