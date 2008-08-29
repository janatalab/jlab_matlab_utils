function status = wt_meg_anal_hdr(pfid,pversion,nfile,all_NChan,NMeg,NReference,NEeg,NAnalog,NBad_chan,all_bad_chan,NEpoch,Epoch,nfreq,obs_labels);

fwrite(pfid,pversion,'integer*2');
fwrite(pfid,nfile,'integer*2');
fwrite(pfid,all_NChan,'integer*2');
fwrite(pfid,NMeg,'integer*2');
fwrite(pfid,NReference,'integer*2');
fwrite(pfid,NEeg,'integer*2');
fwrite(pfid,NAnalog,'integer*2');
fwrite(pfid,NBad_chan,'integer*2');
fwrite(pfid,all_bad_chan,'integer*2');
fwrite(pfid,NEpoch,'integer*2');
fwrite(pfid,Epoch,'real*4');
fwrite(pfid,nfreq,'integer*2');
fwrite(pfid,size(obs_labels,2),'integer*2');
fwrite(pfid,obs_labels,'char*1')
status = 1;
