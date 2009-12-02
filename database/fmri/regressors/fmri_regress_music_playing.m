function [names,vals] = fmri_regress_music_playing(pinfo,minfo,sess)

% generates 'music playing' block-type regressors for fmri data
% 
%   [names,vals] = fmri_regress_music_playing(pinfo,minfo,sess)
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
names = {'music_playing'};
vals = [];

% get presentation data, stimulus onsets
sfilt.include.all = minfo.response_filter;
sinfo = ensemble_filter(pinfo,sfilt);
sc = set_var_col_const(sinfo.vars);
onsets = sinfo.data{sc.RUN_REL_TIME}/1000;

% Durations and amplitude
durs = ones(size(onsets))*minfo.music_dur;
amp  = ones(size(onsets));

% Now build the regressor
vals = fmri_convolve_regress(onsets,durs,amp,...
    pinfo.scanner.TR,pinfo.scanner.dt,pinfo.scanner.actual_nvol);
