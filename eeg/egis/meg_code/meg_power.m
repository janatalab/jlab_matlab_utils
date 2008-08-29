function power_data = meg_power(ravgcsdm,iavgcsdm,NChan,bad_chan,NReference,noise_flag);
mask = ones(1,NChan);
mask(bad_chan) = zeros(1,size(bad_chan,2));
chpair = channel_pairs(NChan);
if strcmp(noise_flag,'no')
	power_data = zeros(size(ravgcsdm,1),NChan);
	for ichan = 1:NChan
			power_data(:,ichan) = ravgcsdm(:,chpair(ichan,ichan));
	end;

	testreal = isreal(sum(power_data));
	if testreal ~= 1
		error('complex numbers in power array')
	end;
else
	power_data = zeros(size(ravgcsdm,1),NChan);
	for ichan = 149:NChan
			power_data(:,ichan) = ravgcsdm(:,chpair(ichan,ichan));
	end;

	testreal = isreal(sum(power_data));
	if testreal ~= 1
		error('complex numbers in power array')
	end;
	for ichan = 1:148
		if mask(ichan)
		chan_list = [ichan 149:151];
		for ifreq = 1:size(ravgcsdm,1)
			for ilist = 1:size(chan_list,2)
				for jlist = 1:size(chan_list,2)
					aug_csdm(ilist,jlist) = ravgcsdm(ifreq,chpair(ilist,jlist))+i*iavgcsdm(ifreq,chpair(ilist,jlist));
					if ilist > jlist
						aug_csdm(ilist,jlist) = conj(aug_csdm(ilist,jlist));
					end;
				end;
			end;
			gan = det(aug_csdm);
			gnn = det(aug_csdm(2:size(aug_csdm,1),2:size(aug_csdm,2)));
			cohan = 1 - (gan/(ravgcsdm(ifreq,chpair(ichan,ichan))*gnn));
			power_data(ifreq,ichan) = (1-cohan)*ravgcsdm(ifreq,chpair(ichan,ichan));
		end;
		end;
	end;
end;

	
	
