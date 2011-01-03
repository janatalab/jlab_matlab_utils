function [names,vals] = fmri_regress_test(pinfo,minfo,sess)

% test function for fmri regressor generation scripts
% 
%   [names,vals] = fmri_regress_test(pinfo,minfo,sess)
% 
% 2009.11.05 FB

names = {'regone','regtwo','regthree'};
vals = rand(50,3);
