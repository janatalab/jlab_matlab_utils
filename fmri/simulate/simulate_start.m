% simulate_start.m
%
% Top-level script for batch processing
% 

% 11/12/03 Petr Janata

if strmatch('MAC',computer)
  warning off MATLAB:warn_parse_feval_script
end

global dataroot batch_root spm_path SCAN_OFFSET sinfo

% List of functions to perform
%
% For each function flag, list subjects for which the function should be run.
% The subject indices refer to indices in the sinfo structure returned by
% get_prime_long_sinfo
%

start_time = clock;

CONVERT_FORMAT = []; % 1:n

RUN_SPM_BATCH = 1;
batch_root = pwd;
batch_file = fullfile(batch_root,'simulate_batch.m');

spm_path = '/usr/local/matlab/toolbox/local/spm99/';

dataroot = '/data1/simulate/';


%
% Get info on all subjects
%

sinfo = get_simulate_sinfo;

%
% Convert data from GE format into AVW
%

if CONVERT_FORMAT
  convert_format(sinfo(CONVERT_FORMAT),dataroot,dataroot, SCAN_OFFSET);
end

%
%  Run the SPM batch file
%

if RUN_SPM_BATCH
  spm_bch(batch_file);
end

elapsed_time = etime(clock,start_time)
