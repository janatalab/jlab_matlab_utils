% get_exp_info.m
%
% See get_exp_info_block.m for more annotations throughout the script
%

% 2003, Petr Janata

ntt = 0;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L1_F2L1';			
tinfo.num_trials(ntt) = 12;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L1_F2L2';			
tinfo.num_trials(ntt) = 12;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L2_F2L1';			
tinfo.num_trials(ntt) = 12;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L2_F2L2';			
tinfo.num_trials(ntt) = 12;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L3_F2L1';			
tinfo.num_trials(ntt) = 6;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L3_F2L2';			
tinfo.num_trials(ntt) = 6;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L4_F2L1';		
tinfo.num_trials(ntt) = 6;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L4_F2L2';		
tinfo.num_trials(ntt) = 6;

ntt = ntt+1;
tinfo.id{ntt} = 'Sil';			
tinfo.num_trials(ntt) = 0; 

trials_per_iteration = sum(tinfo.num_trials);

num_iterations = 2;

nslices = 18;
min_ms_per_slice = 100;

TR = nslices*min_ms_per_slice;

notes_per_seq = 8;
ms_per_note = 700;
resp_period = 2000;
prime_dur = (notes_per_seq-1)*ms_per_note;

cue_dur = 500;
cue_stim_soa_range = [1000 2000];
inter_trial_interval_range = [1*TR 3*TR]*1;

mag_stabilize_s = 30;

total_trials = trials_per_iteration * num_iterations
stim_dur = notes_per_seq*ms_per_note + resp_period;

cue_stim_soa_list = ...
    rand(1,total_trials)*diff(cue_stim_soa_range)+min(cue_stim_soa_range);
trial_durs = cue_stim_soa_list + stim_dur;
pure_stim_time = sum(trial_durs);

iti_list = rand(1,total_trials)*diff(inter_trial_interval_range)+min(inter_trial_interval_range);
total_jitter_time = sum(iti_list);

total_time_s = (pure_stim_time+total_jitter_time)/1000
total_time_min = total_time_s/60

% Construct trial orders for each iteration

trial_type_list = [];
for itt = 1:ntt
  trial_type_list(end+1:end+tinfo.num_trials(itt)) = itt;
end % for itt

% Randomly permute the trial type list
trial_list = trial_type_list(randperm(length(trial_type_list)));

% Calculate onsets
onsets = cumsum(iti_list + trial_durs)-(iti_list(1)+trial_durs(1));

scans_per_iteration = ceil((total_time_s*1000/num_iterations)/TR)