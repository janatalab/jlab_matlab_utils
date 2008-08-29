function [meg_indices, meg_channels, bad_sensors, eeg_indices, eeg_channels,analog_indices, analog_channels, reference_indices, reference_channels] = parse_names(Names, NMeg, NEeg, NBad_Sensors,NAnalog,NReference);
%
%pulls out the index numbers on a trial of data that would be the 
%MEG data and EEG data respectively as meg_indices and eeg_indices
%then creates meg_assign and eeg_assign which have the info to order 
%the data correctly.  Similary there are analog and noise information 
%
%e.g. to get the meg_data ordered correctly 
%
% trialdata(meg_channels) = trialdata(meg_indices);
% trialdata(reference_channels) = trialdata(reference_indices);
% etc. 

trigger_index = find(Names(:,1) == 'T')';
response_index = find(Names(:,1) == 'R')';
utility_index = find(Names(:,1) == 'U')';
meg_indices = find(Names(:,1) == 'A')';
reference_indices = find(Names(:,1) == 'M'|Names(:,1) == 'G')';
bad_sensor_indices = find(Names(:,1) =='E' & Names(:,2) == 'A')';
eeg_indices = find(Names(:,1) == 'E' & Names(:,2) ~= 'A')';
analog_indices = find(Names(:,1) == 'X')';

for i = 1:NMeg
meg_channels(i) = str2num(Names(meg_indices(i),2:5));
end;
total_NMeg = 148;
reference_channels(1:NReference) = [total_NMeg+1:total_NMeg+NReference];

if NEeg > 0
	for i = 1:NEeg
		eeg_channels(i) = str2num(Names(eeg_indices(i),2:5))+103+NReference;
	end;
else
	eeg_channels = [];
end;
if NAnalog > 0
	for i = 1:NAnalog
		analog_channels(i) = str2num(Names(analog_indices(i),2:5))+total_NMeg+NEeg+NReference;
	end;
else
		analog_channels =  [];
end;

if NBad_Sensors > 0
for i = 1:NBad_Sensors
	bad_sensors(i) = str2num(Names(bad_sensor_indices(i),3:6));
end;
else
	bad_sensors = [];
end

if NMeg + NBad_Sensors < 148
	mask = ones(1,148);
	mask(meg_channels) = zeros(1,size(meg_channels,2));
	mask(bad_sensors) = zeros(1,size(bad_sensors,2));
	bad_sensors = [bad_sensors find(mask)];
end;





























