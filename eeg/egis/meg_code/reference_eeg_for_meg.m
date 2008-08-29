function trialdata = reference_eeg_for_meg(trialdata,eeg_channels,mask,EEG_reference);

if strcmp(EEG_reference,'mastoids');

	temp = trialdata(:,eeg_channels) - 0.5*trialdata(:,eeg_channels(4))*ones(1,size(eeg_channels,2));
	trialdata(:,eeg_channels) = temp;
elseif strcmp(EEG_reference,'bipolar1');
%
% bipolars are along one dimension in anterior-posterior directions
% there will be other bipolars for other montages
%
	temp = zeros(size(trialdata,1),size(eeg_channels,2));
	temp(:,1) = trialdata(eeg_channels);
%
% to be finished later
%

end;


	