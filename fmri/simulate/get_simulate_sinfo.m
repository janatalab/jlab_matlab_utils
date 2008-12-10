function s = get_prime_base_sinfo
%
%   sinfo = get_prime_base_sinfo;
%
% NOTE: For all subjects condition orders will be [1 2 3 4]
%

% SCAN_OFFSET:  Specifies the number of images at the beginning of each run
% that should be tossed out.  If event times are specified relative to slice
% acquisition onset, then they must be adjusted by the number of volumes that
% have been jettisoned. This is done in or1_model_protos.m
%
% The volumes are removed in the file conversion stage
%

global SCAN_OFFSET

SCAN_OFFSET = 15; 

s = struct('id','','nexams',[],'exam_nums', [], 'cond_order', [], 'series_mappings',[]);

cond_names = {'Block1', 'Block2', 'Block3', 'Block4'};
vols_per_cond = [191 191 191 191];

nsub = 0;

%
%  Subject 1: 
%

nsub = nsub+1;
s(nsub).id = 'test';
s(nsub).exam_nums = [];
s(nsub).nexams = length(s(nsub).exam_nums);

s(nsub).cond_order = [1 2 3 4]; % if directory 'epi', ifnot 1 2
s(nsub).nvol = vols_per_cond(s(nsub).cond_order);

% enter series mappings for each exam
s(nsub).series_mappings{1} =  [{'002'}, {'coplanar'};
  {'003'} {'epi_12'}];
s(nsub).series_mappings{2} = [{'002'}, {'coplanar2'};
  {'003'} {'epi_34'};
  {'004'} {'hires'}];

