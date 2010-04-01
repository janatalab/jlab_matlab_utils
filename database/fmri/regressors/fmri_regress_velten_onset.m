function [names,vals] = fmri_regress_velten_onset(pinfo,minfo,sess)

% generates fmri data regressors for onset of velten sentence stimuli
% 
%   [names,vals] = fmri_regress_velten_onset(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array of six regressor names
%   vals = volume X regressor matrix containing motion regressors
% 
% FB 2010.03.16

% init output vars
names = {};
vals = [];

if isempty(strfind(pinfo.presfname,'velten')), return, end

onfilt.include.all.EVENT_CODE = {'pic_*'};
onfilt.include.all.EVENT_TYPE = {'Picture'};
ondata = ensemble_filter(pinfo,onfilt);
pc = set_var_col_const(pinfo.vars);
ons = ondata.data{pc.RUN_REL_TIME}/1000;
nons = length(ons);

if ~nons
  fprintf(1,'%s regressor skipped: no onsets found\n',regid);
  return
end

names = {pinfo.regid};

% calculate durations, amplitude
amp = ones(1,nons);
durs = amp';

% Now build the regressor
vals = fmri_convolve_regress(ons,durs,amp,pinfo.scanner.TR,...
    pinfo.scanner.dt,pinfo.scanner.actual_nvol);      
