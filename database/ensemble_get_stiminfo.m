function stimulusDataStruct = ensemble_get_stiminfo(indata,params)
% Obtains stimulus and attribute meta data for stim IDs in the response table.
%
% ensemble_get_stiminfo({respDataStruct,stimids},params)
%
% The first input argument can be either
% respDataStruct - the response table data structure containing stim IDs to retrieve meta data
% stimids - vector of stimulus IDs to obtain information on
% 
%

% 9 Feb, 2007 - First Version, Stefan Tomic
% 15 Mar, 2007 - Edited function to conform to ensemble data structures, S.T.
% 05/18/07 PJ - Added conn_id support
% 8/13/07 ST - minor fix to make vars a row cell array instead of
%              column cell array
% 10/07/08 PJ - enabled passing in of stimids to make this function more useful
% 09/30/09 ST - Added code to deal with missing attribute data (simply set
%               attributeMeta to empty struct)
  
% Handle database and connection defaults

try database = params.ensemble.database; catch database = 'ensemble_main'; end

try
  conn_id = params.ensemble.conn_id;
catch
  conn_id = 1;
  tmp_conn_id = 1;
end

if isstruct(indata)
  %obtain a list of unique stimulus IDs that are in the response table
  stimIDCol = strmatch('stimulus_id',indata.vars,'exact');
  stimIDCol = unique( stimIDCol(~isnan(stimIDCol)) );

  %throw an error if no stimulus IDs were provided
  if(isempty(stimIDCol))
    error(['Cannot get stimulus information. Stimulus IDs not provided' ...
	  ' from response data']);
  end

  stimIDAll = indata.data{stimIDCol};
  stimIDList = unique(stimIDAll(~isnan(stimIDAll)));
else
  stimIDList = unique(indata);
end

% Extract the stimulus data
stimMeta = mysql_extract_metadata('table','stimulus', ...
    'stimulus_id',stimIDList, ...
    'conn_id', conn_id);

%add the stimulus_id field to attribute structures, so that we can
%separate the attribute meta info from the stimulus meta info
for iStim = 1:length(stimMeta)
    for iAttribute = 1:length(stimMeta(iStim).attribute)
      stimMeta(iStim).attribute(iAttribute).stimulus_id = stimMeta(iStim).stimulus_id;
    end
end

%separate the stimulus and attribute meta info
attributeMeta = [stimMeta.attribute];
if(~isempty(attributeMeta))
  attributeMeta = orderfields(attributeMeta,{'stimulus_id','attribute_id','name','class','attribute_value_double','attribute_value_text'});
else
  %if there is no attribute metadata, just setting this to an empty struct
  attributeMeta = struct('stimulus_id',[],'attribute_id',[],'name',[],'class',[],'attribute_value_double',[],'attribute_value_text',[]);
end
stimMeta = rmfield(stimMeta,'attribute');

att2DCells = squeeze(struct2cell(attributeMeta))';
stim2DCells = squeeze(struct2cell(stimMeta))';

%reorganize the 2D cells to cell column data, and set type and vars
%fields
stimMetaCells = ensemble_init_data_struct;
stimMetaCells.name = 'stimulus_metadata';
stimMetaCells.type = 'stimulus_metadata';
for iStimCol = 1:size(stim2DCells,2)
  if isnumeric(stim2DCells{1,iStimCol})
    columnData = [stim2DCells{:,iStimCol}]';
  else
    columnData = {stim2DCells{:,iStimCol}}';
  end
  stimMetaCells.data(:,iStimCol) = {columnData};
end
stimMetaCells.vars = fieldnames(stimMeta)';

attMetaCells = ensemble_init_data_struct;
attMetaCells.name = 'stimulus_x_attribute_metadata';
attMetaCells.type = 'stimulus_x_attribute_metadata';
for iAttCol = 1:size(att2DCells,2)  
  if isnumeric(att2DCells{1,iAttCol})
    columnData = [att2DCells{:,iAttCol}]';
  else
    columnData = {att2DCells{:,iAttCol}}';
  end
  attMetaCells.data(:,iAttCol) = {columnData};
end
attMetaCells.vars = fieldnames(attributeMeta)';

stimulusDataStruct = ensemble_init_data_struct;
stimulusDataStruct.vars = {'stimulus_metadata', ...
		    'stimulus_x_attribute_metadata'};
stimulusDataStruct.data = {stimMetaCells, attMetaCells};

if(exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
