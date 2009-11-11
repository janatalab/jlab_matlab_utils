function [names,vals] = fmri_regress_stim(pinfo,minfo,sess)

% generates fmri data regressors for stimulus-based responses
% 
%   [names,vals] = fmri_regress_stim(pinfo,minfo,sess)
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

regid = pinfo.regid;

% init output vars
names = {regid};
vals = [];

% init vars
ccm = minfo.cond_cue_map; % function handle for condition to cue mappings
cue = ccm(regid); % get cue for this regressor
resp_params = extract_resp_params_v2(cue,pinfo,minfo,sess);
pc = set_var_col_const(pinfo.vars);

if isempty(resp_params)
  return
end

% Now get the vector of stimulus onset times to which these response
% parameters refer
sfilt.include.all = minfo.response_filter;
sinfo = ensemble_filter(pinfo,sfilt);
onsets = sinfo.data{pc.RUN_REL_TIME}/1000;

% Durations
durs = ones(size(onsets))*minfo.music_dur;

% Now build the regressor
vals = fmri_convolve_regress(onsets,durs,resp_params,pinfo.scanner.TR,...
    pinfo.scanner.dt,pinfo.scanner.actual_nvol);