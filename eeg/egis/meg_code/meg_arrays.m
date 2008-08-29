function [arrays_meg] = meg_arrays(array_label,bad_chan);
[xelec, yelec, zelec] = sensors(148);
if strcmp(array_label,'left')
     arrays_temp = find(yelec > 0);
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
elseif strcmp(array_label,'right')
     arrays_temp = find(yelec < 0);
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
elseif strcmp(array_label,'posterior')
     arrays_temp = find(xelec < 0);
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
elseif strcmp(array_label,'anterior')
     arrays_temp = find(xelec >= 0);
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
elseif strcmp(array_label,'left_lat')
     arrays_temp = [113 95 139 121 103 138 120 102 137 119 101 136 118 100 135 117 99 134 133 116 98 132 115 97 131 114 96];
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
elseif strcmp(array_label,'right_lat')
     arrays_temp = [140 122 104 141 123 105 142 124 106 143 125 144 126 108 145 127 109 146 128 110 147 129 111 148 130 112];
     temp = zeros(1,148);
     temp(arrays_temp) = ones(1,length(arrays_temp));
     temp(bad_chan) = zeros(1,length(bad_chan));
     arrays_meg = find(temp);
else
     error('unknown array');
end;

