% get_exp_info.m
%
% Specifies some basic properties of an fMRI experiment design matrix. The default
% version of this script accommodates a slightly mutated 2 factor design with
% multiple levels in each factor. However, the number of factors/levels that enter
% into the design matrix can easily be modified either at this stage or in simulate_model_protos.m
%
% This script was originally written to prepare for a music cognition experiment
% in which there were trials with primes and targets.  Some of the nomenclature
% for the calculations of prime durations has been retained.  Simply think of
% 'notes' as 'events' in a sequence of priming events preceding a target.
%
% This script will also spit back estimates of the experiment duration.
%

% 2003, Petr Janata

clear tinfo

ntt = 0;

% Note: In this example, Factor 1 refers to the prime type and Factor 2 refers
% to an attribute of the target.  Other target attributes exist in relation to
% the prime, but these are handled at the next level of model construction, not
% at the level of trial type specification.  
%
% If you are not interested in prime/target trial structure, it is easy to
% eliminate one or the other element by forcing the durations of the respective
% element to zero and then not including that element as a condition in the
% model (in simulate_model_protos.m).

ntt = ntt+1;
tinfo.id{ntt} = 'F1L1_F2L1';  % Factor 1, Level 1; Factor 2 Level 1
tinfo.num_trials(ntt) = 2;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L1_F2L2';  % Factor 1, Level 1; Factor 2 Level 2	
tinfo.num_trials(ntt) = 2;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L2_F2L1';  % Factor 1, Level 2; Factor 2 Level 1
tinfo.num_trials(ntt) = 2;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L2_F2L2';	% Factor 1, Level 2; Factor 2 Level 2	
tinfo.num_trials(ntt) = 2;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L3_F2L1';	% Factor 1, Level 3; Factor 2 Level 1	
tinfo.num_trials(ntt) = 0;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L3_F2L2';	% Factor 1, Level 3; Factor 2 Level 2	
tinfo.num_trials(ntt) = 0;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L4_F2L1';	% Factor 1, Level 4; Factor 2 Level 1	
tinfo.num_trials(ntt) = 0;

ntt = ntt+1;
tinfo.id{ntt} = 'F1L4_F2L2';	% Factor 1, Level 4; Factor 2 Level 2	
tinfo.num_trials(ntt) = 0;

ntt = ntt+1;
tinfo.id{ntt} = 'Sil';			% Silence (blank trials)
tinfo.num_trials(ntt) = 0; 

trials_per_iteration = sum(tinfo.num_trials);

num_iterations = 1;  % Number of iterations through the number of trials
                     % specified above

nslices = 18;				% number of slices to acquire
min_ms_per_slice = 100;  % The fastest sampling rate we can go at, i.e. 10 slices/s

TR = nslices*min_ms_per_slice;		% The experiment's TR

notes_per_seq = 1;  % Number of events in priming sequence
ms_per_note = 30000;  % Duration (ms) of each event in priming sequence
resp_period = 15000*2;  % Amount of time (ms) allowed for response to target
prime_dur = (notes_per_seq-1)*ms_per_note;  % Duration of prime

cue_dur = 500;  % Duration of trial type cue
cue_stim_soa_range = [1000 2000];	% Range of asynchronies (ms) between
                                        % trial type cue and trial onset
inter_trial_interval_range = [1*TR 3*TR]*0; % range of times between trials

total_trials = trials_per_iteration * num_iterations
stim_dur = notes_per_seq*ms_per_note + resp_period; % overall trial duration

% Create randomized list of trial type cue to stimulus onset asynchronies
cue_stim_soa_list = ...
    rand(1,total_trials)*diff(cue_stim_soa_range)+min(cue_stim_soa_range);
trial_durs = cue_stim_soa_list + stim_dur;
pure_stim_time = sum(trial_durs);

% Create randomized list of inter-trial intervals
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

% Calculate trial onsets
onsets = cumsum(iti_list + trial_durs)-(iti_list(1)+trial_durs(1));

% Number of volumes collected for each iteration through the trials
scans_per_iteration = ceil((total_time_s*1000/num_iterations)/TR)
