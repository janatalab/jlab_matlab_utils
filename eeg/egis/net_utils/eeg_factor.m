function [factor_mag, factor_ph, factor_eig] = eeg_factor(avgcsdm,good_chan,freq,output);
ch_pair_indices	
if size(avgcsdm,2) == 8385
	NChan = 129;
else
	error('number of channels is unknown');
end;

if any(good_chan) > 129
	error('channel out of range');
end;

if freq > size(avgcsdm,1)
	error('frequency out of range');
end;
csdm_data = avgcsdm;
csdm_mat = zeros(NChan,NChan);
csdm = zeros(size(good_chan,2),size(good_chan,2));
for ichan =1:NChan
	for jchan = ichan:NChan;
		csdm_mat(ichan,jchan) = csdm_data(freq,chpair(ichan,jchan));
		csdm_mat(jchan,ichan) = conj(csdm_mat(ichan,jchan));
	end;	
end;
	
csdm = csdm_mat(good_chan,good_chan);

[vec_csd, val_csd] = eig(csdm,'nobalance');
val_csd = real(val_csd);
factor_eig = diag(val_csd);
vec_csd = vec_csd';
factor_mag = zeros(size(vec_csd,1),129);
factor_mag(:,good_chan) = abs(vec_csd);
phase_data = vec_csd;
if output == 'cos'
     phase_data = real(phase_data)./abs(phase_data);
else
                                real_sign = sign(real(phase_data));
                                phase_data = angle(phase_data);
                                ang_sign = sign(phase_data);
                                for iph = 1:size(phase_data,1)
                                        for jph = size(phase_data,2)
                                                if real_sign(iph,jph) < 0
                                                        phase_data(iph,jph) = phase_data(iph,jph) + pi*real_sign(iph,jph)*ang_sign(iph,jph);
                                                end;
                                        end;
                                end;


                              
end
			phase_data(find(isnan(phase_data))) = zeros(size(find(isnan(phase_data))));

			if ~isreal(phase_data)
                  		error('Phase matrix is not real')
                	end

factor_ph = zeros(size(phase_data,1),129);
factor_ph(:,good_chan) = phase_data;



