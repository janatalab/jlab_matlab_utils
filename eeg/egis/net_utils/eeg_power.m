function power_data = eeg_power(ravgcsdm)
%power_data = eeg_power(ravgcsdm)
%
%extracts power values from avgcsdm matrix
%
ch_pair_indices;
if size(ravgcsdm,2) == 8385
	NChan = 129;
else
	error('unknown number of channels')
end;
power_data = zeros(size(ravgcsdm,1),NChan);
for ichan = 1:NChan
	power_data(:,ichan) = ravgcsdm(:,chpair(ichan,ichan));
end;

testreal = isreal(sum(power_data));

if testreal ~= 1
	error('complex numbers in power array')
end;


