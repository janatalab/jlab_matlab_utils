function EEG = rename_events(EEG, ei);
% EEG = rename_events(EEG, event_info);
%
% Matches event codes in 'code' field of event_info with event codes in the 'type' field of the
% EEG.event structure (see EEGLAB for more info about the EEG structure), and
% replaces the value in the type field with the 'name' field in event_info
%

% August 7, 2005 Petr Janata - initiated script

nevents = length(ei);

for ie = 1:nevents
  % Get a list of indices that match this event
  % How we do this depends on the variable types in event.type
  
  if any(cellfun('isclass',{EEG.event.type},'char'))
    event_idxs = find(ismember({EEG.event.type},num2str(ei(ie).code)));
  else
    event_idxs = find(ismember([EEG.event.type],ei(ie).code));
  end
  nidx = length(event_idxs);
  
  if nidx
    % Get the latencies and modify them so that they turn out right at the end of
    % the process
    tmp = num2cell([EEG.event(event_idxs).latency]/EEG.srate);
  
    % Create the newdata array
    clear newdata
    [newdata{1:nidx,1}] = deal(ei(ie).name);
    [newdata{:,2}] = deal(tmp{:});
  
    assignin('base','newdata',newdata);
    
    % Make the call to the EEGLAB function that will make the modificationss
    warning off
    EEG = pop_importevent(EEG,'event','newdata','indices',event_idxs, ...
	'fields',{'type','latency'},'append','yes');
    warning on
  end % if nidx
end % for ie

return
