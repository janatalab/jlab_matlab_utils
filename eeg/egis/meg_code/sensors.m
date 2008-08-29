function [xsensor,ysensor,zsensor] = sensors(NMeg);

if NMeg ~= 148
	error('improper number of MEG sensors');
end;

load '/home/ramesh/matlab_code/meg_code/sensor_loc.asc';

xsensor = sensor_loc(:,1)';
ysensor = sensor_loc(:,2)';
zsensor = sensor_loc(:,3)';

