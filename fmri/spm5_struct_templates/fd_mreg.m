function [jobs] = fd_mreg
% Specification of a multiple regression design
% The tree below is incomplete (FB)
% -jobs 
%  \-jobs{1}
%    \-stats 
%      \-stats{1}
%        \-factorial_design
%          |-des
%          | \-mreg
%          |   \-scans ............ {}
%          |   \-incint ........... 1  (include the intercept?)
%          |   \-mcov
%          |     \-c .............. []
%          |     \-cname .......... ''
%          |     \-iCC ............ 1  (1: center on overall mean, 0:don't)
%          | 
%          |-masking
%          | |-tm
%          | | \-tm_none ........... []
%          | |-im ................ 1-x-1 double
%          | \-em 
%          |   \-em{1} ............. []
%          |-globalc
%          | \-g_omit ............ []
%          |-globalm
%          | |-gmsca
%          | | \-gmsca_no .......... []
%          | \-glonorm ........... 1-x-1 double
%          \-dir ............... <UNDEFINED>
% -jobs 
%    \-stats 
% 
% FB <fbarrett@ucdavis.edu> 2012.02.29 - based on fd_ff.m

% Factorial design
jobs{1}.stats{1}.factorial_design.des.mreg = struct(...
    'scans','<UNDEFINED>', ...
    'incint',1, ...
    'mcov',struct('c',[],...
            'cname','',...
            'iCC',1));

% Covariates
jobs{1}.stats{1}.factorial_design.cov = struct([]);

% Masking
jobs{1}.stats{1}.factorial_design.masking.tm.tm_none = ...
	reshape(double([]),[0,0]);
jobs{1}.stats{1}.factorial_design.masking.im = reshape(double(1),[1,1]);
jobs{1}.stats{1}.factorial_design.masking.em{1} = char('');

jobs{1}.stats{1}.factorial_design.globalc.g_omit = reshape(double([]),[0,0]);
jobs{1}.stats{1}.factorial_design.globalm.gmsca.gmsca_no = reshape(double([]),[0,0]);
jobs{1}.stats{1}.factorial_design.globalm.glonorm = reshape(double(1),[1,1]);
jobs{1}.stats{1}.factorial_design.dir = char('<UNDEFINED>');
