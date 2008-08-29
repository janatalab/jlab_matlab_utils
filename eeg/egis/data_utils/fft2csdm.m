function cross_trial = fft2csdm(fft_trialdata)

cross_trial = zeros(fix(size(fft_trialdata,1)),(size(fft_trialdata,2).^2+size(fft_trialdata,2))/2); 
icount = 1;

if strcmp(computer, 'MAC2')
  for i = 1:size(fft_trialdata,2)
	for j = i:size(fft_trialdata,2)
		cross_trial(:,icount) = ppc_cmplx_mult(fft_trial(1:fix(size(fft_trialdata,1)),i),conj(fft_trial(1:size(fft_trialdata,1),j)));
		icount = icount+1;
	end;
  end;
else
  for i = 1:size(fft_trialdata,2)
	for j = i:size(fft_trialdata,2)
		cross_trial(:,icount) = fft_trialdata(1:fix(size(fft_trialdata,1)),i).*conj(fft_trialdata(1:fix(size(fft_trialdata,1)),j));
		icount = icount+1;
	end;
  end;
end

