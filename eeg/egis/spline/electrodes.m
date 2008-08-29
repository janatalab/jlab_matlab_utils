function [xelec,yelec,zelec] = electrodes(NChan)
global code_path
global position_file
position_file = ['''' code_path  'oregon129.xyz'''];
if NChan == 129
eval(['load ' position_file ' -ascii']);
oregon = oregon129';
xelec = oregon(2,:);
yelec = oregon(3,:);
zelec = oregon(4,:);
end;

