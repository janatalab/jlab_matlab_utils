function [version,NChan_file,NMeg,NReference,NEeg,NAnalog,NBad_chan,bad_chan,NEpoch,Epoch,nfreq] = rd_meg_csdm_hdr(csdmfname);
%
% reads just the header
%

if isempty(csdmfname)
	[fid,csdmfname] = get_fid('rb');
	fclose(fid);
end;

csdmfid = open_file_w_byte_order(csdmfname,-1);

version = fread(csdmfid,1,'integer*2');
NChan_file = fread(csdmfid,1,'integer*2');
NMeg = fread(csdmfid,1,'integer*2');
NReference = fread(csdmfid,1,'integer*2');
NEeg = fread(csdmfid,1,'integer*2');
NAnalog = fread(csdmfid,1,'integer*2');
NBad_chan = fread(csdmfid,1,'integer*2');
bad_chan = fread(csdmfid,[1,NBad_chan],'integer*2');
NEpoch = fread(csdmfid,1,'integer*2');
Epoch = fread(csdmfid,1,'real*4');
nfreq = fread(csdmfid,1,'integer*2');
fclose(csdmfid);


