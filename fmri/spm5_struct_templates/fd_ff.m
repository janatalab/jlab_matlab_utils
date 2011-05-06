function [jobs] = fd_ff
% Specification of full-factorial design
% The tree below is incomplete (PJ)
% -jobs 
%  \-jobs{1}
%    \-stats 
%      \-stats{1}
%        \-factorial_design
%          |-des
%          | \-fd
%              \-fact
%          |   \-icell
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

% Factorial design
jobs{1}.stats{1}.factorial_design.des.fd.fact = struct(...
	'name', '<UNDEFINED>', ...
	'levels', '<UNDEFINED>', ...
	'dept', 1, ...  % 1=repeated-measures/non-sphericity
	'variance', 1, ...
	'gmsca', 0, ...
	'ancova', 0);

jobs{1}.stats{1}.factorial_design.des.fd.icell = struct(...
	'levels', '<UNDEFINED>', ...
	'scans', '<UNDEFINED>');
		 
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
