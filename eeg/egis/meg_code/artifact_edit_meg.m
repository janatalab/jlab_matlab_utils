function mask = artifact_edit_meg(trialdata,meg_channels,eeg_channels,MEG_max,EEG_max,add_bad_chan);
%
% applies amplitude criterion to MEG and EEG channels to artifact edit. 
% returning a mask with ones for good channels and 0 for bad channels 
%

mask = ones(1,size(trialdata,2));

max_trialdata = max(abs(trialdata));
bad_meg = find(max_trialdata(meg_channels) > MEG_max);
if ~isempty(eeg_channels)
	bad_eeg = find(max_trialdata(eeg_channels) > EEG_max);
	mask(eeg_channels(bad_eeg)) =ones(size(bad_eeg,2));
end;

mask(meg_channels(bad_meg)) = zeros(1,size(bad_meg,2));
if ~isempty(add_bad_chan)
	mask(add_bad_chan) = zeros(1,size(add_bad_chan,2));
end;

