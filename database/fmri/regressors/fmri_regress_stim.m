function [names,vals] = fmri_regress_stim(pinfo,minfo,sess)

% generates fmri data regressors for stimulus-based responses
% 
%   [names,vals] = fmri_regress_stim(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% also handles limits based on timespan judgments
% NOTE: need to add a little of 
% 
% REQUIRES (incomplete list)
%   pinfo
%   minfo
%       .stim
%           .(regid)
%   sess
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
ccm = parse_fh(minfo.cond_cue_map); % function handle for condition to cue mappings
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
if isfield(minfo,'music_dur')
  % Durations are constant
  durs = ones(size(onsets))*minfo.music_dur;
elseif isfield(minfo,'music_dur_db') && minfo.music_dur_db
  sids = sinfo.data{pc.EVENT_CODE};
  durs = fmri_stim_duration(pinfo,minfo,sids);
end

% limit by timespan?
if ~isempty(strfind(regid,'timespan'))
  [onsets,durs,resp_params] = ...
      fmri_regress_timespan(pinfo,minfo,onsets,durs,resp_params);
end % if ~isempty(strfind(regid,'timespan

% special instructions for this regid?
if isfield(minfo,'stim') && isfield(minfo.stim,regid)
  rinfo = minfo.stim.(regid);
  
  % generate dummy regressors for different levels of this stim?
  if isfield(rinfo,'dummy_regs')
    %%% THIS CONDITIONAL RETURNS, IF SET, NOTHING AFTER THIS WILL BE RUN 
    sv = rinfo.dummy_regs.seg_vals;
    ns = length(sv);
    if isfield(rinfo.dummy_regs,'seg_names')
      sn = rinfo.dummy_regs.seg_names;
    else
      sn = {};
      for k=1:ns, sn{k} = sprintf('%s_%d',regid,k); end
    end
    vals = [];
    names = {};
    for k=1:length(sv)
      li = ismember(resp_params,sv{k});
      if isempty(li), continue, end
      lval = fmri_convolve_regress(onsets(li),durs(li),ones(length(li),1),...
          pinfo.scanner.TR,pinfo.scanner.dt,pinfo.scanner.actual_nvol);
      vals = [vals lval];
      names = [names sn{k}];
    end % for k=1:length(sv)
    
    %%%% return!!!
    return
  end % if isfield(rinfo,'dummy_regs
end % if isfield(minfo,'stim') && isfield(minfo.stim,regid

% Now build the regressor
vals = fmri_convolve_regress(onsets,durs,resp_params,pinfo.scanner.TR,...
    pinfo.scanner.dt,pinfo.scanner.actual_nvol);
