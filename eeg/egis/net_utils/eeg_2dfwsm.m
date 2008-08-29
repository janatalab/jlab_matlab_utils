function [fwsm_mat] = eeg_2dfwsm(csdm_data,good_chan,freq,k_range);

ch_pair_indices

fwsm_mat = zeros(size(k_range,2),size(k_range,2));
NChan = 129;
%for if = 1:size(freqs,2)

csdm_mat = zeros(129,129);
csdm = zeros(size(good_chan,2),size(good_chan,2));
for ichan =1:NChan
	for jchan = ichan:NChan;
		csdm_mat(ichan,jchan) = csdm_data(freq,chpair(ichan,jchan));
		csdm_mat(jchan,ichan) = conj(csdm_mat(ichan,jchan));
	end;	
end;

csdm = csdm_mat(good_chan,good_chan);
[x,y,z] = electrodes(129);
elec = xyz2tp(x,y,z);
cchan_pos(:,1) = 9.2*elec(:,1);
cchan_pos(:,2) = 9.2*elec(:,2);
chan_pos = cchan_pos(good_chan,:);
for c = 1:size(k_range,2)
for j = 1:size(k_range,2)
	kvec1 = exp(i*k_range(j)*(ones(size(good_chan,2),1)*chan_pos(:,2)'-chan_pos(:,2)*ones(1,size(good_chan,2))));
	kvec2 = exp(i*k_range(c)*(ones(size(good_chan,2),1)*chan_pos(:,1)'-chan_pos(:,1)*ones(1,size(good_chan,2))));
	kvec = kvec1.*kvec2;
%	for k = 1:size(good_chan,2)
%		for l = 1:size(good_chan,2) 
%			kvec(k,l) = exp(i*k_range(j)*(chan_pos(l,2)-chan_pos(k,2)))*exp(i*k_range(c)*(chan_pos(l,1)-chan_pos(k,1)));
%					
%	end;
%	end;

	power = csdm.*kvec;
	fwsm_mat(c,j) = real(sum(sum(power))/(size(good_chan,2).^2));
end;
end;







