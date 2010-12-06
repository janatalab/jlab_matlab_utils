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
% FB 2010.02.17 - added function to check databae for stimulus durations

% get presentation data, stimulus onsets
sfilt.include.all = minfo.condlist.music_playing;
sinfo = ensemble_filter(pinfo,sfilt);
sc = set_var_col_const(sinfo.vars);
ons = sinfo.data{sc.RUN_REL_TIME}/1000;
sids = sinfo.data{sc.EVENT_CODE};

if isfield(minfo,'music_dur')
  % Durations are constant
  durs = ones(size(ons))*minfo.music_dur;
elseif isfield(minfo, 'condparams') && ...
		isfield(minfo.condparams, 'music_playing') && ...
		isfield(minfo.condparams.music_playing, 'music_dur') ...
		
	durs = ones(size(ons))*minfo.condparams.music_playing.music_dur;
	
elseif isfield(minfo,'music_dur_db') && minfo.music_dur_db
  durs = fmri_stim_duration(pinfo,minfo,sids);
end
