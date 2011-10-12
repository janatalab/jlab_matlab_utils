function sessData = ensemble_session_trial_info_without_trial_id(indata,params)
% sessData = ensemble_session_trial_info_without_trial_id(indata,params)
% 
% Extracts the response_order and stimulus IDs for each given session.
%
% Adapted from ensemble_session_trial_info to accomodate data from
% experiments for which trial IDs were not created
%
% indata is a cell array of two data structures (response_data and
% session_info), as returned from ensemble_load_expinfo
% 
% This function looks at the response table (from ensemble_load_expinfo)
% of an experiment to extract the response_order and stimulus IDs for each 
% session, in the order that they were presented. This info is stored 
% in a 'trial info' structure that is tagged to the end of each session 
% in the session structure (also retrieved from ensemble_load_expinfo). 
%
% Optionally,the audio corresponding to each trial can be parsed from a session recording
% (e.g. from digital performer). The time offsets of each trial relative to
% the beginning of the recording are then recorded. This information is then
% used to parse a continuously recorded MIDI time course into individual responses
% corresponding to each trial (if MIDI data was recorded along-side the audio 
% recording.)
% 
%
% handles only single stimulus trials for now
%
%
% main differences from ensemble_session_trial_info:
%   - accomodates data with no trial_id values in response table
%   - calls parse_midi_slider_response, which resamples midi responses into time courses of evenly spaced points
%   - allows removal of practice trials not recorded in audio but recorded in response_table 
%  
%
%
% 29 June 2011 - BH

if ~isfield(params,'verbose'), params.verbose=0; end

% params.numPracticeTrials parameter is used for removing practice trials from the
% parsed audio/ensemble response set. This is useful if practice
% trials were recorded, but are not desired in the results. Note that this
% should be used only if practice trials recorded in audio file.
if(isfield(params,'numPracticeTrials')) %fix to handle scenario where practice trials not recorded in DP but exist in response table
  practiceTrialIdxs = 1:params.numPracticeTrials;
else
  practiceTrialIdxs = [];    
end

% params.num_practice_trials_database_only paramter is used to remove 
% practice trials from ensemble response set only, and not audio. This is an
% alternative to numPractice Trials and is usefule if practice trials were not
% recorded in the audio files but exist in the response table
if(isfield(params,'num_practice_trials_database_only'))
    databasePracticeTrialIdxs = 1:params.num_practice_trials_database_only;
else
    databasePracticeTrialIdxs = [];
end

dataStructCrit.name = 'session_info';
sessIdx = ensemble_find_analysis_struct(indata,dataStructCrit);
sessData = indata{sessIdx};
sessData.name = 'session_x_trial_info';
sessData.type = 'session_x_trial_info';
sessData.vars{end+1} = 'trial_info';
sessCols = set_var_col_const(sessData.vars);

dataStructCrit.name = 'response_data';
respIdx = ensemble_find_analysis_struct(indata,dataStructCrit);
respData = indata{respIdx};
respCols = set_var_col_const(respData.vars);


sessionIDs = sessData.data{sessCols.session_id};
subjectIDs = sessData.data{sessCols.subject_id};

for sessIdx = 1:length(sessionIDs)

  sessInfo.session_id = sessionIDs(sessIdx);
  sessInfo.subject_id  = subjectIDs{sessIdx};
 
  message(['Sorting trial info for session # ' num2str(sessIdx) ...
	   ', ID ' sessInfo.subject_id ...
       ', Ensemble session # ' num2str(sessInfo.session_id)],...
       params.verbose);
  
  filt.include.all.session_id = sessInfo.session_id;
  respThisSess = ensemble_filter(respData,filt);
  
  trialInfoStruct = ensemble_init_data_struct;
  trialInfoStruct.name='trial info';
  trialInfoStruct.type='trial info';
  trialInfoStruct.vars={'response_order' 'stimulus_id'};

   
  % get rid of repeated response order values and stimulus ids
  % this is done by taking the diff of the response order
  % and zeroing out all the indices where the diff == 0
  % also, getting rid of NaNs in the and stim id list 
  trialInfoVector = [respThisSess.data{respCols.response_order}  ...
      respThisSess.data{respCols.stimulus_id}];
  
  trialdiffIdx = [ 1 diff(trialInfoVector(:,1))'];
  trialInfoVector(trialdiffIdx == 0,:) = [];
  trialInfoVector( isnan(trialInfoVector(:,2)),: ) = [];
  
  % if databasePracticeTrialIdxs set remove those trials
  if ~isempty(databasePracticeTrialIdxs)
      trialInfoVector(databasePracticeTrialIdxs,:) = [];
  end
  
  trialInfoCols = set_var_col_const(trialInfoStruct.vars);
  trialInfoStruct.data{trialInfoCols.response_order}=trialInfoVector(:,1);
  trialInfoStruct.data{trialInfoCols.stimulus_id}=trialInfoVector(:,2);
  
  %if params.parse_audio was set, then audio parsing will be
  %performed. params.parse_audio serves as the param struct for the
  %audio parser
  if isfield(params,'parse_audio_stims') && ~isempty(params.parse_audio_stims)
  
      parseAudioParams = params.parse_audio_stims;
      parseAudioParams.filename = ...
          replaceFilenameTags(parseAudioParams.filename,sessInfo);


      %if the recorded responses in Ensemble or audio recordings
      %don't match, these parameters will ignore either Ensemble
      %session info (ignoreMatchedEventForSub) or audio info
      %(ignoreParsedStimAudioForSub)during matching
      if(isfield(params,'ignoreMatchedEventForSub'))
          
          [hasEventExclusions,ignoreCellIdx] = ismember(sessInfo.subject_id,params.ignoreMatchedEventForSub{1});
          if(ignoreCellIdx ~= 0)
              parseAudioParams.ignoreEventIdxs = ...
                  params.ignoreMatchedEventForSub{2}{ignoreCellIdx};


              %if ensemble practice sessions are missing from database,
              %then we might need to revise which practice trials to
              %throw out.
              practiceTrialIdxs = setdiff(practiceTrialIdxs,parseAudioParams.ignoreEventIdxs);

          end
      end
        
      if(isfield(params,'ignoreParsedStimAudioForSub'))
          [hasStimAudioExclusions,ignoreCellIdx] = ...
          ismember(sessInfo.subject_id,params.ignoreParsedStimAudioForSub{1});
          
          if(ignoreCellIdx ~= 0)
              parseAudioParams.ignoreParsedAudioIdxs = params.ignoreParsedStimAudioForSub{2}{ignoreCellIdx};
          end
      end
      
      trialInfoStruct = parse_audio_clips(trialInfoStruct,parseAudioParams);
      if isempty(trialInfoStruct)
          warning('no audio data for %s',sessInfo.subject_id);
          continue
      end
  end %parse audio
  
  if isfield(params,'parse_midi_resps') && ~isempty(params.parse_midi_resps)
    
    parseMidiRespsParams = params.parse_midi_resps;
    parseMidiRespsParams.filename = replaceFilenameTags(parseMidiRespsParams.filename,sessInfo);
    trialInfoStruct = parse_midi_slider_response(trialInfoStruct,parseMidiRespsParams);
   
  end %parse MIDI
  
  %throw out the practice trials if this was set
  allTrialIdxs = 1:length(trialInfoStruct.data{1});
  trialIdxsToKeep = setdiff(allTrialIdxs,practiceTrialIdxs);
  for iCol = 1:length(trialInfoStruct.data)
    trialInfoStruct.data{iCol} = trialInfoStruct.data{iCol}(trialIdxsToKeep);
  end
  
  sessData.data{sessCols.trial_info}(sessIdx,1) = {trialInfoStruct};
  
end %sessIdx

function message(messageString,verbose)

if(verbose > 0)
  disp(sprintf(messageString));
end

return


function filename = replaceFilenameTags(filename,replaceStrings)

filename = strrep(filename,'<subject_id>',replaceStrings.subject_id);
filename = strrep(filename,'<session_id>',num2str(replaceStrings.session_id));

return