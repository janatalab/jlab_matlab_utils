function ordered_trial = fix_trial_order(trialdata,meg_channels,meg_indices,eeg_channels,eeg_indices,analog_channels,analog_indices,reference_channels,reference_indices);
%
%reorders data into the correct order coming off MEG convert
%
%

ordered_trial(:,meg_channels) = trialdata(:,meg_indices)*(10^12);
ordered_trial(:,eeg_channels) = trialdata(:,eeg_indices)*100;
ordered_trial(:,analog_channels) = trialdata(:,analog_indices)*100;
ordered_trial(:,reference_channels) = trialdata(:,reference_indices)*(10^12);

