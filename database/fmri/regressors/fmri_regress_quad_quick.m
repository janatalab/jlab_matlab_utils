function [names,vals] = fmri_regress_quad_quick(pinfo,minfo,sess)

% generates a quickly increasing quadratic run regressor for fmri data
% 
%   [names,vals] = fmri_regress_quad_quick(pinfo,minfo,sess)
% 
% called by fmri_generate_regress
% 
% REQUIRES
% 
% RETURNS
%   names = cell array containing one regressor name: 'linear_run%d'
%   vals = -(log(t/T)/log(t(1)/T)) + 1
% 
% FB 2009.11.05

names = {sprintf('quad_slow_run%d',pinfo.irun)};
vals = -(log((1:pinfo.scanner.actual_nvol)./pinfo.scanner.actual_nvol)/...
    log(1/pinfo.scanner.actual_nvol)) + 1;
