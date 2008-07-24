% or1_start.m
%
% Top-level script for analyzing data from Reg's object categorization experiment
%
% 11/22/00 Petr Janata
% 02/12/01 PJ -- implemented volume tossing during conversion step

global dataroot spm_path SCAN_OFFSET

% List of functions to perform
%
% For each function flag, list subjects for which the function should be run.
% The subject indices refer to indices in the sinfo structure returned by
% get_ts1_sinfo
%

start_time = clock;

CONVERT_FORMAT = [2]; %

RUN_SPM_BATCH = 0;
batch_file = '/data1/matlab/or1/or1_batch.m';

diary_name = sprintf('or1_diary_%s_%s', datestr(datenum(now),1), datestr(datenum(now),15));
%diary(diary_name);

NOTIFY = 0;

%indata_root = '/mnt/cdrom/';
indata_root = '/data1/anatomicals';
outdata_root = '/data1/anatomicals';
spm_root = '/data1/or1/';

dataroot = spm_root;

spm_path = '/usr/local/matlab/toolbox/local/spm99/';

%
% Get info on all subjects
%

s = get_anat_sinfo;

%
% Convert data from GE format into AVW
%

if CONVERT_FORMAT
  convert_format(s(CONVERT_FORMAT),indata_root,outdata_root, SCAN_OFFSET);
end

%
%  Run the SPM batch file
%

if RUN_SPM_BATCH
  spm_bch(batch_file);
end

%diary on
elapsed_time = etime(clock,start_time)

if NOTIFY
  unix(['mail janata@dartmouth.edu -s "or1_batch finished" < /dev/null'])
%  unix('mail babsi@dartmouth.edu -s "ts1_batch finished" < /dev/null')
end

