function [onsets,durs] = fmri_regress_cues(pinfo,minfo,sess)

% generates cue onset regressors
% 
%   [names,vals] = fmri_regress_cues(pinfo,minfo,sess)
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
onsets=[];
durs=[];

cfilt.include.all = minfo.condlist.cues;
lcond = ensemble_filter(pinfo,cfilt);
lc = set_var_col_const(lcond.vars);

onsets = lcond.data{lc.RUN_REL_TIME}/1000;
durs = zeros(size(onsets));
