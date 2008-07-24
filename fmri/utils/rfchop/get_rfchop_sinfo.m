function s = get_rfchop_sinfo
% get_rfchop_sinfo

global dataroot 

dataroot = '/data1/rfchop/';

s = struct('id','','nexams',[],'exam_nums', [], 'cond_order', [], 'series_mappings',[]);

cond_names = {'Phantom_RFon', 'Phantom_RFoff', 'Finger_Run1', 'Finger_Run2', 'Finger_Run3', 'Finger_Run4'};
vols_per_cond = [150 150 86 86 86 86];

nsub = 0;

nsub = nsub+1;
s(nsub).id = 'RFchoptest';
s(nsub).sub_id = 1;
s(nsub).exam_nums = [2122];
s(nsub).nexams = length(s(nsub).exam_nums);

s(nsub).cond_order = [1 2];
s(nsub).use_runs = [1 2];
s(nsub).nvol = vols_per_cond(s(nsub).cond_order);

s(nsub).series_mappings{1} =  [{'002'}, {'epi'}, {1};
  {'003'} {'epi'} {2}];

nsub = nsub+1;
s(nsub).id = '08may02PJ';
s(nsub).sub_id = 2;
s(nsub).exam_nums = [2123];
s(nsub).nexams = length(s(nsub).exam_nums);

s(nsub).cond_order = [3:6];
s(nsub).use_runs = [1:4];
s(nsub).nvol = vols_per_cond(s(nsub).cond_order);

s(nsub).series_mappings{1} =  [{'002'}, {'coplanar'}, {[]};
  {'003'}, {'epi'}, {1};
  {'004'} {'epi'} {2};
  {'005'} {'epi'} {3};
  {'006'} {'epi'} {4}];
