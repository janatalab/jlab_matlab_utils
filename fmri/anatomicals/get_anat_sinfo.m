function s = get_anat_sinfo
%
%   sinfo = get_anat_sinfo;
%
% Information on various subjects for whom anatomical images were acquired
% outside the context of an fMRI experiment, e.g. for EEG coregistration
%
					
s = struct('id','','nexams',[],'exam_nums', [], 'cond_order', [], 'series_mappings',[]);

cond_names = {};

nsub = 0;

%
%  Subject 1: 
%

vols_per_cond = [];

nsub = nsub+1;
s(nsub).id = '19jun01JB';
s(nsub).exam_nums = [1200];
s(nsub).nexams = length(s(nsub).exam_nums);

s(nsub).cond_order = [];
s(nsub).nvol = vols_per_cond(s(nsub).cond_order);

% enter series mappings for each exam
s(nsub).series_mappings{1} = ...
    [{'003'}, {'hires'}]

%
%  Subject 2: 
%

vols_per_cond = [];

nsub = nsub+1;
s(nsub).id = '19jun01BT';
s(nsub).exam_nums = [1201];
s(nsub).nexams = length(s(nsub).exam_nums);

s(nsub).cond_order = [];
s(nsub).nvol = vols_per_cond(s(nsub).cond_order);

% enter series mappings for each exam
s(nsub).series_mappings{1} = ...
    [{'002'}, {'hires'}]

