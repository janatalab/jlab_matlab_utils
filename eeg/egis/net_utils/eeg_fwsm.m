function [fwsm_mat] = eeg_fwsm(csdm_data,good_chan,chan_pos,freqs,k_range);

  ch_pair_indices

  fwsm_mat = zeros(size(freqs,2),size(k_range,2));
  NChan = 129;
  for if = 1:size(freqs,2)

    csdm_mat = zeros(129,129);
    csdm = zeros(size(good_chan,2),size(good_chan,2));
    for ichan =1:NChan
      for jchan = ichan:NChan;
	csdm_mat(ichan,jchan) = csdm_data(freqs(if),chpair(ichan,jchan));
	csdm_mat(jchan,ichan) = conj(csdm_mat(ichan,jchan));
      end;	
    end;

    csdm = csdm_mat(good_chan,good_chan);

    for j = 1:size(k_range,2)
      for k = 1:size(good_chan,2)
	for l = 1:size(good_chan,2) 
	  kvec(k,l) = exp(i*k_range(j)*(chan_pos(l)-chan_pos(k)));
	  
	end;
      end;

      power = csdm.*kvec;
      fwsm_mat(if,j) = real(sum(sum(power))/(size(good_chan,2).^2));
    end;
  end;
