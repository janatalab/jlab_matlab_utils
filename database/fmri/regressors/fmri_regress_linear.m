function [names,vals] = fmri_regress_linear(pinfo,minfo,sess)

% generates a linear run regressor for fmri data
% 
%   [names,vals] = fmri_regress_linear(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array containing one regressor name: 'linear_run%d'
%   vals = nvolX1 monotonically increasing regressor, ranging from -1 to 1    
% 
% FB 2009.11.05

names = {sprintf('linear_run%d',pinfo.irun)};
vals = [0:1/(pinfo.scanner.actual_nvol-1):1]';
