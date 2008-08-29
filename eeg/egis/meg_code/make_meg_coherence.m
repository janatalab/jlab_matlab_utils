function status = make_meg_coherence(csdmfname,coherencefname,obs_labels,noise_flag);
%status = make_meg_coherence(csdmfname,coherencefname,obs_labels,noise_reduction);
%
%csdmfname = matrix of csdmfilenames
%obs_labels = labels for the csdm matrices
%noise_flag = 'yes' or 'no'
%	 if yes the reference channels are used, and the mutiple 
%	coherence between reference channels and each sensor is 
%	calculated.  This tells us the proportion of the coherence spectrum
%	due to the noise. We then remove this. See Bandat and Piersol, 1986
%

if nargin < 3
	error('not enough input arguments');
end;

if nargin == 3
	noise_flag = 'yes';
end;

if isempty(obs_labels)
	obs_labels = csdmfname;
end;
all_bad_chan = zeros(size(csdmfname,1),100);
for i = 1:size(csdmfname,1)
	[version, NChan_file,NMeg(i),NReference(i),NEeg(i),NAnalog(i),NBad_chan(i),bad_chan,NEpoch(i),Epoch(i),nfreq(i)] = rd_meg_csdm_hdr(csdmfname(i,:));
	all_bad_chan(i,1:NBad_chan(i)) = bad_chan;
	all_NChan(i) = 148+NReference(i)+NEeg(i) +NAnalog(i);
end;

pfid = fopen(coherencefname,'wb');
pversion = -3;
nfile = size(csdmfname,1);
status = wt_meg_anal_hdr(pfid,pversion,nfile,all_NChan,NMeg,NReference,NEeg,NAnalog,NBad_chan,all_bad_chan,NEpoch,Epoch,nfreq,obs_labels);


for i = 1:nfile
	[ravgcsdm,iavgcsdm] = rd_meg_csdm(csdmfname(i,:));
	coherence_data = meg_coherence(ravgcsdm,iavgcsdm,all_NChan(i),all_bad_chan(i,1:NBad_chan(i)),NReference(i),noise_flag);
	fwrite(pfid,coherence_data','real*4');
end;

fclose('all')
status = 1;
	




