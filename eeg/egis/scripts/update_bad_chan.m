function status = update_bad_chan(new_bad_chan,powerfname);

if nargin < 2
	[fid,powerfname] = get_fid('rb');
else
	fid = fopen(powerfname,'rb');
end;
outfid = fopen([powerfname '.fix'],'wb');
version = fread(fid,1,'int16')
fwrite(outfid,version,'int16');
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
obs_labels = setstr(obs_labels);
sbad = size(new_bad_chan,2);
for i = 1:size(bad_chan,1)
	bad_chan(i,Nbad_chan(i)+1:Nbad_chan(i)+sbad) = new_bad_chan;
	Nbad_chan(i) = Nbad_chan(i)+sbad;
end;

wt_anal_hdr(outfid,nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq);

for i = 1:nfiles
	power = fread(fid,[NChan(i), NFreq(i)],'real*4');
	fwrite(outfid,power,'real*4');
end;
fclose('all')
status = 1;	


