function [names,vals] = fmri_regress_average(pinfo,minfo,sess)

% generates an averaged fmri regressor from a given set of responses
% 
%   [names,vals] = fmri_regress_average(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% note: expects the same number of responses for all variables in a given
% set
% 
% REQUIRES
%   minfo.stim_average.var_sets - cell array of cell array of strings
%       containing cues to be averaged together
%   minfo.stim_average.set_names - cell array of strings that are names for
%       each of the sets defined above
% 
% RETURNS
%   names = cell array of six regressor names
%   vals = volume X regressor matrix containing motion regressors
% 
% FB 2010.02.25

% init output vars
names = minfo.stim_average.set_names;
nsets = length(names);
vals = zeros(pinfo.scanner.TR,nsets);
stim_avg_vars = minfo.stim_average.var_sets;
nsav = length(stim_avg_vars);
if nsets ~= nsav, error('different number of sets and names'); end

% init vars
ccm = parse_fh(minfo.cond_cue_map); % function handle for condition to cue mappings

% Now get the vector of stimulus onset times to which these response
% parameters refer
sfilt.include.all = minfo.response_filter;
sinfo = ensemble_filter(pinfo,sfilt);
pc = set_var_col_const(pinfo.vars);
onsets = sinfo.data{pc.RUN_REL_TIME}/1000;

% Durations
if isfield(minfo,'music_dur')
  % Durations are constant
  durs = ones(size(onsets))*minfo.music_dur;
elseif isfield(minfo,'music_dur_db') && minfo.music_dur_db
  sids = sinfo.data{pc.EVENT_CODE};
  durs = fmri_stim_duration(pinfo,minfo,sids);
end

% iterate over stim avg sets, get average resp_params, calculate vals
vals = zeros(pinfo.scanner.actual_nvol,nsav);
for j=1:nsav
  vset = stim_avg_vars{j};
  resp_params = [];
  for k=1:length(vset)
    cue = ccm(vset{k}); % get cue for this regressor
    lrp = extract_resp_params_v2(cue,pinfo,minfo,sess);
    if isempty(lrp), continue; end
    resp_params = [resp_params lrp];
  end % for k=1:length(vset
  avg_rp = mean(resp_params,2);
  vals(:,j) = fmri_convolve_regress(onsets,durs,avg_rp,pinfo.scanner.TR,...
      pinfo.scanner.dt,pinfo.scanner.actual_nvol);
end % for j=1:length(stim_avg_vars
