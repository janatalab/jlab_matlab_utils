function status = wt_anal_hdr(pfid,nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,Nfreq);

fwrite(pfid,nfiles,'int16');
fwrite(pfid,size(obs_labels,2),'int16');
for i = 1:size(obs_labels,1)
	fwrite(pfid,obs_labels(i,:),'char');
end;
fwrite(pfid,Epoch,'int16');
fwrite(pfid,Window_Length,'int16')
fwrite(pfid,NEpoch,'int16');	
fwrite(pfid,Nbad_chan,'int16');
fwrite(pfid,bad_chan,'int16');
fwrite(pfid,ref_flag,'int16');
fwrite(pfid,reference,'char');		
fwrite(pfid,NChan,'int16');
fwrite(pfid,Nfreq,'int16');
status = 1;