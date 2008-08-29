function phase_data = eeg_phase(ravgcsdm, iavgcsdm,output)
%phase_data= eeg_phase(ravgcsdm,iavgcsdm,output); 
%
%calculates relative phase 
%
%avgcsdm= cross spectral density matrix
%rel_chan = channel for indexing relative phase
%output = output type can be 'pha' or 'cos' for raw phase and cosine of phase
%  defaults to 'cos' (optional)
%
if nargin < 2
	error('duuh')
end;
if nargin < 3
	output = 'cos';
end;
ch_pair_indices;
phase_data = zeros(size(ravgcsdm,1),size(ravgcsdm,2));
if strcmp(output,'cos')
	phase_data = ravgcsdm./sqrt(ravgcsdm.^2 + iavgcsdm.^2);
else
%	real_sign = sign(ravgcsdm);
	phase_data = atan2(iavgcsdm,ravgcsdm);
%	ang_sign = sign(phase_data);
%	for iph = 1:size(phase_data,1)
%		for jph = size(phase_data,2)
%			if real_sign(iph,jph) < 0
%				phase_data(iph,jph) = phase_data(iph,jph) + pi*real_sign(iph,jph)*ang_sign(iph,jph);
%			end;
%		end;
%	end;
%			
end
phase_data(find(isnan(phase_data))) = zeros(size(find(isnan(phase_data))));
testreal = isreal(sum(phase_data));
if testreal ~= 1
	error('complex numbers')
end;


