function [names,vals] = fmri_regress_motion_per_run(pinfo,minfo,sess)

% generates separate motion regressors for each run of fmri data
% 
%   [names,vals] = fmri_regress_motion_per_run(pinfo,minfo,sess)
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
names = cell(1,6);
vals = zeros(pinfo.scanner.actual_nvol,6);

% init vars
motvars = {'Motion-X','Motion-Y','Motion-Z','Motion-Roll','Motion-Pitch','Motion-Yaw'};

% load motion params
X = load(pinfo.motionparam_fname);

% iterate over dimensions, extract data
for ivar=1:6
  names{ivar} = sprintf('%s_%d',motvars{ivar},pinfo.irun);
  tmp = X(1:pinfo.scanner.actual_nvol,ivar);
  vals(:,ivar) = tmp-mean(tmp); % remove mean
end
