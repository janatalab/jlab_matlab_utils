function [cross_trial,fft_trial] = csdm(trialdata)

trialdata = zeromean(trialdata);
fft_trial = fft(trialdata)./(size(trialdata,1)+1);
cross_trial = zeros(fix(size(trialdata,1)/2),(size(trialdata,2).^2+size(trialdata,2))/2); 
icount = 1;

if strcmp(computer, 'MAC2')
  for i = 1:size(trialdata,2)
	for j = i:size(trialdata,2)
		cross_trial(:,icount) = ppc_cmplx_mult(fft_trial(1:fix(size(trialdata,1)/2),i),conj(fft_trial(1:size(trialdata,1)/2,j)));
		icount = icount+1;
	end;
  end;
else
  for i = 1:size(trialdata,2)
	for j = i:size(trialdata,2)
		cross_trial(:,icount) = fft_trial(1:fix(size(trialdata,1)/2),i).*conj(fft_trial(1:fix(size(trialdata,1)/2),j));
		icount = icount+1;
	end;
  end;
end
