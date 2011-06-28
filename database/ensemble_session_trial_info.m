function sessData = ensemble_session_trial_info(indata,params)
% Extracts the trial and session IDs in order, for each given session.
%
% indata is a cell array of two data structures (response_data and
% session_info), as returned from ensemble_load_expinfo
%
% This function looks at the response table (from ensemble_load_expinfo)
% of an experiment to extract the trial IDs and stimulus IDs for each 
% session, in the order which they were presented. This info is stored 
% in a 'trial info' structure that is tagged to the end of each session 
% in the session structure (also retrieved from ensemble_load_expinfo). Optionally,
% the audio corresponding to each trial can be parsed from a session recording
% (e.g. from digital performer). The time offsets of each trial relative to
% the beginning of the recording are then recorded. This information is then
% used to parse MIDI responses corresponding to each trial (if MIDI data was
% recorded along-side the audio recording.
% 
%
% handles only single stimulus trials for now
%
%
% 21 March 2007 - S.T.

% params.numPracticeTrials parameter is used for removing practice trials from the
% parsed audio/ensemble response set. This is useful if practice
% trials were recorded, but are not desired in the results
if(isfield(params,'numPracticeTrials'))
  practiceTrialIdxs = 1:params.numPracticeTrials;
else
  practiceTrialIdxs = [];
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
	   ', ID ' num2str(sessInfo.session_id)],params.verbose);
 
  
  filt.include.all.session_id = sessInfo.session_id;
  respThisSess = ensemble_filter(respData,filt);

  trialInfoStruct = ensemble_init_data_struct;
  trialInfoStruct.name='trial info';
  trialInfoStruct.type='trial info';
  trialInfoStruct.vars={'trial_id' 'stimulus_id'};
  
  % get rid of repeated trial ids and stimulus ids
  % this is done by taking the diff of the response order
  % and zeroing out all the indices where the diff == 0
  % also, getting rid of NaNs in the trial id and stim id list
  trialInfoVector = [respThisSess.data{respCols.response_order}  ...
      respThisSess.data{respCols.trial_id}  ...
      respThisSess.data{respCols.stimulus_id}];
  
  trialdiffIdx = [ 1 diff(trialInfoVector(:,1))'];
  trialInfoVector(trialdiffIdx == 0,:) = [];
  trialInfoVector( isnan(trialInfoVector(:,2)),: ) = [];
  
  trialInfoCols = set_var_col_const(trialInfoStruct.vars);
  trialInfoStruct.data{trialInfoCols.trial_id}=trialInfoVector(:,2);
  trialInfoStruct.data{trialInfoCols.stimulus_id}=trialInfoVector(:,3);
 
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
  
  end
  
  if isfield(params,'parse_midi_resps') && ~isempty(params.parse_midi_resps)
    
    parseMidiRespsParams = params.parse_midi_resps;
    parseMidiRespsParams.filename = replaceFilenameTags(parseMidiRespsParams.filename,sessInfo);
    trialInfoStruct = parse_midi_responses(trialInfoStruct,parseMidiRespsParams);
   
  end
  
  %throw out the practice trials if this was set
  allTrialIdxs = 1:length(trialInfoStruct.data{1});
  trialIdxsToKeep = setdiff(allTrialIdxs,practiceTrialIdxs);
  for iCol = 1:length(trialInfoStruct.data)
    trialInfoStruct.data{iCol} = trialInfoStruct.data{iCol}(trialIdxsToKeep);
  end
  
  sessData.data{sessCols.trial_info}(sessIdx,1) = {trialInfoStruct};
  
end

function message(messageString,verbose)

if(verbose > 0)
  disp(sprintf(messageString));
end

return


function filename = replaceFilenameTags(filename,replaceStrings)

filename = strrep(filename,'<subject_id>',replaceStrings.subject_id);
filename = strrep(filename,'<session_id>',num2str(replaceStrings.session_id));

return
