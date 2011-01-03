function [names,vals] = fmri_regress_quad_slow(pinfo,minfo,sess)

% generates a slowly increasing quadratic run regressor for fmri data
% 
%   [names,vals] = fmri_regress_quad_slow(pinfo,minfo,sess)
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

names = {sprintf('quad_slow_run%d',pinfo.irun)};
vals = [((1:pinfo.scanner.actual_nvol)/pinfo.scanner.actual_nvol).^2]';
