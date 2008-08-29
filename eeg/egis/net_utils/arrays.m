function [channel_arrays,group_names] = arrays(array_flag,bad_chan);
%[channel_arrays,group_names] = arrays(array_flag,bad_chan)
%
% electrode arrays for the 129 channel EEG system
% argin
% array_flag = 'oned','quadrant','hemisphere','all','frontback'
% bad_chan = bad channel list
%
% argout
% channel_arrays = arrays of channels named array_flag. Rows are distinct lists
%	the arrays are zero padded to have the same lenghth. 
% group_names = names for arrays
if strcmp(array_flag,'oned')
	[xelec,yelec,zelec] = electrodes(129);
	twod_dist = twod_pos(xelec,yelec,zelec,9.2);
	ch_pair_indices;
	channel_arrays = [18 23 27 34 40 46 50 58 65 71;19 24 28 35 41 47 51 59 66 72;20 25 29 36 42 48 52 60 67 0; 12 21 30 37 43 53 61 0 0 0; 16 11 6 129 55 62 68 73 0 0; 5 119 112 105 94 87 79 0 0 0; 4 124 118 111 104 99 93 86 78 0; 10 3 123 117 110 103 98 92 85 77; 15 9 2 122 116 109 102 97 91 84];
	group_names = ['L1';'L2';'L3';'L4'; 'M0';'R4';'R3';'R2';'R1'];
elseif strcmp(array_flag,'quadrant')
	
	channel_arrays = [1 2 3 4 5 8 9 10 14 15 121 122 123 124 119 113 107 116 117 118 110 111 112 104 105 106;81 88 94 99 103 80 87 93 98 102 79 86 92 97 101 78 85 91 96 77 84 90 83 0 0 0; 32 38 43 48 47 46 50 51 52 53 54 57 58 59 60 61 64 65 66 67 70 71 72 75 0 0; 22 18 23 19 26 27 24 20 12 33 39 34 28 25 21 7 35 29 30 31 36 37 40 41 42 0];
	group_names = ['RA'; 'RP'; 'LP'; 'LA'];

elseif strcmp(array_flag,'hemisphere')
	channel_arrays = [1 2 3 4 5 8 9 10 14 15 121 122 123 124 119 113 107 116 117 118 110 111 112 104 105 106 81 88 94 99 103 80 87 93 98 102 79 86 92 97 101 78 85 91 96 77 84 90 83; 32 38 43 48 47 46 50 51 52 53 54 57 58 59 60 61 64 65 66 67 70 71 72 75 22 18 23 19 26 27 24 20 12 33 39 34 28 25 21 7 35 29 30 31 36 37 40 41 42];
	group_names = ['R'; 'L'];
	
elseif strcmp(array_flag, 'all')
	channel_arrays = [1 2 3 4 5 8 9 10 14 15 121 122 123 124 119 113 107 116 117 118 110 111 112 104 105 106 81 88 94 99 103 80 87 93 98 102 79 86 92 97 101 78 85 91 96 77 84 90 83 32 38 43 48 47 46 50 51 52 53 54 57 58 59 60 61 64 65 66 67 70 71 72 75 22 18 23 19 26 27 24 20 12 33 39 34 28 25 21 7 35 29 30 31 36 37 40 41 42];
	group_names = 'A';
elseif strcmp(array_flag, 'frontback')
	channel_arrays2 = [1 2 3 4 5 8 9 10 14 15 121 122 123 124 119 113 107 116 117 118 110 111 112 104 105 106;81 88 94 99 103 80 87 93 98 102 79 86 92 97 101 78 85 91 96 77 84 90 83 0 0 0; 32 38 43 48 47 46 50 51 52 53 54 57 58 59 60 61 64 65 66 67 70 71 72 75 0 0; 22 18 23 19 26 27 24 20 12 33 39 34 28 25 21 7 35 29 30 31 36 37 40 41 42 0];
	group_names = ['AN'; 'PO'];
	channel_arrays = [channel_arrays2(1,:) channel_arrays2(4,:); channel_arrays2(2,:) channel_arrays2(3,:)];
elseif strcmp(array_flag, 'aug_oned')
	channel_arrays = [19 24 28 35 41 47 51 59 66 72 20 25 29 36 42 48 52 60 67 0 0 0; 12 21 30 37 43 53 61 16 11 6 129 55 62 68 73 5 119 112 105 94 87 79; 4 124 118 111 104 99 93 86 78 10 3 123 117 110 103 98 92 85 77 0 0 0];
	group_names = ['L';'M';'R'];
else 
	error('unknown array');

end;

for i = 1:size(channel_arrays,1)
	temp = zeros(1,size(find(channel_arrays(i,:)),2));
	temp = channel_arrays(i,find(channel_arrays(i,:)));
	mask = zeros(1,129);
	mask(temp) = ones(1,size(temp,2));
	mask(bad_chan) = zeros(1,size(bad_chan,2));
	temp2 = zeros(1,size(find(mask),2));
	temp2 = find(mask);
	channel_arrays(i,:) = zeros(1,size(channel_arrays,2));
	channel_arrays(i,1:size(temp2,2)) = temp2;
end;



