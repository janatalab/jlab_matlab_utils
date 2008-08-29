function [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,Nfreq] = rd_anal_hdr(fid);

nfiles = fread(fid,1,'int16');
length_obs = fread(fid,1,'int16');
for i = 1:nfiles
	obs_labels(i,:) = fread(fid,[1,length_obs],'char');
end;
Epoch = fread(fid,[1,nfiles],'int16');
Window_Length = fread(fid,[1,nfiles],'int16');
NEpoch = fread(fid,[1,nfiles],'int16');
Nbad_chan = fread(fid,[1,nfiles],'int16');
bad_chan = fread(fid,[nfiles,129],'int16');
ref_flag = fread(fid,[1,nfiles],'int16');
reference = fread(fid,[nfiles,129],'char');
NChan = fread(fid,[1,nfiles],'int16');
Nfreq = fread(fid,[1,nfiles],'int16');

