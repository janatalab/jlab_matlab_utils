% simulate_model_protos.m
%
% Creates a design matrix
%
% In combination with the get_exp_info* scripts, a number of different types of
% experiment designs can be simulated with relative ease, ranging from block
% designs to event designs, including cue-target type of designs.  Originally,
% the script was written to model a priming experiment using musical material,
% but that's just semantics. 
%
% In combination with run_post_process.m, this set of scripts can be used to
% examine properties of the design matrix and the repercussions that various
% design choices may have on the evaluation of data.

% 11/12/03 Petr Janata - adapted for simulating various design matrices

global model_action

clear tinfo

nmodels = 1;
nruns = 1;				% this will have to be set dynamically

% Since we may want to try various types of models but not have to revise the
% files containing the trial specifications each time we want to go back and
% forth, we can refer to different models by different names and have the script
% be smart about constructing the design matrix.  Here we create an array of
% structure with experiment information. We can easily add to the array as we add
% experiments we want to simulate.  We then select one of these "experiments"

nexp = 0;

% Create the first element in the exp_info structure array

% increment the experiment list counter
nexp = nexp+1;				

% specify an arbitrary name for the model
exp_info(nexp).name = 'block_1fact1level'; 

% specify the member of the get_exp_info_* family of scripts we are going to call
exp_info(nexp).script_stub = 'block'; % i.e., get_exp_info_block.m

% List the conditions that will be in the model. Note, these names might refer
% to the factors in the model, or factor/level combinations. Another portion of
% this script will associate 'trial types' specified in the get_exp_info.m
% script with each of the conditions in this list
exp_info(nexp).cond_names = {'Tonic'};  % List of conditions we want to evaluate

% Create the next element in the exp_info structure array
nexp = nexp+1;
exp_info(nexp).name = 'block_1fact2levels';
exp_info(nexp).script_stub = 'block';
exp_info(nexp).cond_names = {'Tonic','Subdominant'};  

nexp = nexp+1;
exp_info(nexp).name = 'block_2fact';
exp_info(nexp).script_stub = 'block';
exp_info(nexp).cond_names = {'Tonic','Subdominant','Consonant','Dissonant'};

nexp = nexp+1;
exp_info(nexp).name = 'block_2fact_bylevel';
exp_info(nexp).script_stub = 'block';
exp_info(nexp).cond_names = {'TC','TD','SC','SD'};

nexp = nexp+1;
exp_info(nexp).name = 'event_byfact';
exp_info(nexp).script_stub = 'event';
exp_info(nexp).cond_names = {'Cue silence','Cue music','Tonal','Atonal','Tonic','Subdominant','Random','Consonant','Dissonant'};

nexp = nexp+1;
exp_info(nexp).name = 'event_bylevel';
exp_info(nexp).script_stub = 'event';
exp_info(nexp).cond_names = {'Cue silence','Cue music','Tonal','Atonal','TC','TD','SC','SD','RC','RD'};

nexp = nexp+1;
exp_info(nexp).name = 'event_bylevel_w_silence';
exp_info(nexp).script_stub = 'event';
exp_info(nexp).cond_names = {'Cue silence','Cue music','Silence','Tonal','Atonal','TC','TD','SC','SD','RC','RD'};

nexp = nexp+1;
exp_info(nexp).name = 'event_bylevel_nosil'; % make sure to set number of
                                             % silent trials=0 in get_exp_info_
exp_info(nexp).script_stub = 'event';
exp_info(nexp).cond_names = {'Cue music','Tonal','Atonal','TC','TD','SC','SD','RC','RD'};

% select which of the experiments we just specified is going to be modeled.
use_exp = nexp;  

% Run the script with all of the trial type and timing parameters, etc. for this experiment
script_str = sprintf('get_exp_info_%s', exp_info(use_exp).script_stub);
eval(script_str);

% Build the appropriate condition list
cond_names = exp_info(use_exp).cond_names;
nc = length(cond_names);

% Specify onset times and other info for each condition
% Model everything as events with variable-length durations

bf_ep_idx = zeros(1,nc);
bf_ev_idx = ones(1,nc);

%
% This is the part of the code that specifies which trial types will be associated
% with which conditions (members of the cond_names field in the exp_info
% structure created above) in the design matrix.  All of the timing information will
% get pulled from the structure and variables created in the get_exp_info
% script. Note that some of the conditions model activity associated with the
% prime part of each trial, whereas others reflect activity associated with the
% targets.  The timing info is set accordingly.
%
% You can easily code novel trial_type combinations by adding a case statement
% below.

for ic = 1:nc
  switch cond_names{ic}
    %
    % Conditions that model the trial type cues
    %
    case 'Cue silence'
      % Find all of the desired trial types in the list of trials 
      idxs = find(ismember(trial_list,find(ismember(tinfo.id,'Sil'))));

      % Get the event onset times for all of these trial types
      sot{ic} = onsets(idxs)/TR; 

      % Specify the events are events or epochs. Note: events with durations
      % are effectively epochs
      cond_types{ic} = 'events';
      
      % Specify the event durations in TR units
      durs{ic} = ones(size(idxs))*cue_dur/TR;  

    case 'Cue music'
      idxs = find(ismember(trial_list,find(~ismember(tinfo.id,'Sil'))));
      sot{ic} = onsets(idxs)/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*cue_dur/TR;  % duration in msec
    
    %
    % Conditions that model the type of prime
    %
    case 'Tonal'  % Factor 1, Levels 1 and 2, irrespective of Factor 2 level
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L1','F1L1_F2L2','F1L2_F2L1','F1L2_F2L2'})))); 
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs))/TR; % onset of prime
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*prime_dur/TR; % duration of prime

    case 'Atonal'  % Factor 1, Levels 3 and 4, irrespective of Factor 2 level
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L3_F2L1','F1L3_F2L2','F1L4_F2L1','F1L4_F2L2'})))); 
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs))/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*prime_dur/TR;

    case 'Silence'  % This could be thought of as Factor 1, Level 5
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'Sil'})))); 
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs))/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*prime_dur/TR;     
	
      %
      % Conditions that model the targets
      %
	  
    case 'Tonic'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L1','F1L1_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*ms_per_note/TR;
    
    case 'TC' % tonic, consonant
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L1'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;
	
    case 'TD' % tonic, dissonant
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;

    case 'Subdominant'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L2_F2L1','F1L2_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;

    case 'SC' % subdominant, consonant
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L2_F2L1'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;

    case 'SD' % subdominant, dissonant
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L2_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*ms_per_note/TR;

    case 'Random'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L3_F2L1','F1L3_F2L2','F1L4_F2L1','F1L4_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*ms_per_note/TR;
	  
    case 'RC'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L3_F2L1','F1L4_F2L1'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
      durs{ic} = ones(size(idxs))*ms_per_note/TR;
	
    case 'RD'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L3_F2L2','F1L4_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;
	  
    case 'Consonant'  % By Factor 2, Level 1, irrespective of prime
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L1','F1L2_F2L1','F1L3_F2L1','F1L4_F2L1'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;

    case 'Dissonant'
      idxs = find(ismember(trial_list, ...
	  find(ismember(tinfo.id,{'F1L1_F2L2','F1L2_F2L2','F1L3_F2L2','F1L4_F2L2'}))));
      sot{ic} = (onsets(idxs)+cue_stim_soa_list(idxs) + prime_dur)/TR;
      cond_types{ic} = 'events';
	  durs{ic} = ones(size(idxs))*ms_per_note/TR;

  end % switch
end % for ic=

%
% Fill in the generic model info
%

model_action_list = char({'specify','review','estimate','spec_and_estimate'});

for imodel = 1:nmodels
  model(imodel) = struct( ...
      'types', strmatch(model_action,model_action_list), ... % 1=specify, 2=review, 3=estimate, 4=specify+estimate
      'nsess',          nruns, ...
      'RT',             TR/1000, ...
      'nscans',         [scans_per_iteration], ...	% this is actually specified in the calling script
      'replicated',     0, ...
      'same_time_param', 0, ...
      'conditions_nb',  ones(1,nruns) * nc(imodel), ...     
      'conditions',     1:nc, ...	
      'stochastics_flag', zeros(1,nruns), ...
      'stochastics',    [], ...
      'parametrics_type', {{'none','none'}}, ...
      'parametrics',    [], ...
      'regressors_nb',  zeros(1,nruns), ...
      'regressors',     zeros(1,nruns), ... 
      'global_effects', {'scaling'}, ...
      'burst_mode',     0, ...
      'HF_fil',         'none',  ...
      'HF_cut',         [], ...
      'LF_fil',         'none', ...
      'LF_cut',         [], ...
      'int_corr',       'none', ... 
      'trial_fcon',     0, ...
      'now_later',      1, ...		% 1 = now, 0 = later
      'stop_writing',   0, ...
      'files',          [] ...
      );

%  model_prototype(imodel).nscans = model_prototype(imodel).nscans - SCAN_OFFSET;
end

%
%  Conditions structures
%

for irun = 1:nruns
  conditions(irun) = struct( ...
      'volterra', 0, ...
      'variable_dur', 1, ...
      'names',   {cond_names}, ...
      'onsets',  {sot}, ...   
      'types',    {cond_types}, ... 
      'durations', {durs}, ...
      'bf_ev',   bf_ev_idx, ...
      'bf_ep',   bf_ep_idx ...
      );
end

%
%  Basis function structures
%

bf_ev(1) = struct( ...
  'ev_type', 1, ...			% 1=hrf
  'win_len', [], ...
  'order', [], ...
  'gamma_idxs', [] ...
);

