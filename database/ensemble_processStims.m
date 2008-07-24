function outData = ensemble_processStims(inData,params)
% Wrapper script that applies a given function to a given set of stims.
%
% This is a general wrapper script that cycles through all 
% the stims that are passed in the
% stimulus_metadata structure. Processes to be run on the stims are
% passed in as strings that specify function handles.
% parameters such as stimulus_id,artist,song_name,location are
% automatically passed to the child function.
%
% July 24, 2007 - Stefan Tomic


%placing inData into a cell. This facilitates the following
%analysis search, which currently is not necessary, but will be
%useful if this function later supports other input data structures
%(e.g. stimulus_x_attribute_metadata)
if(~iscell(inData))
  inData = {inData};
end

clear searchCrit;
searchCrit.type = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimMetaData = inData{findIdx};
stimMetaDataCols = set_var_col_const(stimMetaData.vars);

if(isfield(params,'filt'))
  stimMetaData = ensemble_filter(stimMetaData,params.filt);
end

nStims = length(stimMetaData.data{stimMetaDataCols.stimulus_id});

if(~isfield(params,'startStim'))
  params.startStim = 1;
end

for iStim = params.startStim:nStims

  funcParams = params.funcParams;
  
  thisStimID = ...
      stimMetaData.data{stimMetaDataCols.stimulus_id}(iStim);
  funcParams.stimulus_id= thisStimID;
  fprintf('Processing Stim %d, ID %d\n',iStim,thisStimID);
  
  funcParams.artist = stimMetaData.data{stimMetaDataCols.artist}{iStim};
  funcParams.song_name = stimMetaData.data{stimMetaDataCols.name}{iStim};
  funcParams.location = stimMetaData.data{stimMetaDataCols.location}{iStim};
  
  funcH = str2func(params.funcName);
  
  outData{iStim} = funcH(inData,funcParams);
    
end
