function [ons,durs] = fmri_cond_music_playing(pinfo,minfo,sess)

% generates 'music playing' block-type condition definition for fmri data
% 
%   [ons,durs] = fmri_cond_music_playing(pinfo,minfo,sess)
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

% get presentation data, stimulus onsets
sfilt.include.all = minfo.response_filter;
sinfo = ensemble_filter(pinfo,sfilt);
sc = set_var_col_const(sinfo.vars);
ons = sinfo.data{sc.RUN_REL_TIME}/1000;

% Durations and amplitude
durs = ones(size(ons))*minfo.music_dur;
