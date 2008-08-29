function [Epoch,NEpoch,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask] = read_mask(rawfname);
%
%[Epoch,NEpoch,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask] = read_mask(rawfname);
%
%reads in editcode information fom '*.mask' file
%
% rawfname must be a NSFragger filename.  [rawfname '.mask'] must be the mask file, stored in the same directory as the datafile.   
%
maskfname = [rawfname '.mask'];
infid = fopen(maskfname,'rb');
if infid < 0
	error('invalid edit code file name')
end;
Epoch = fread(infid,1,'integer*2');
NEpoch = fread(infid,1,'integer*2');
NChan = fread(infid,1,'integer*2');
Max_microv = fread(infid,1,'integer*2');
Min_Chan = fread(infid,1,'integer*2');
Nbad_chan = fread(infid,1,'integer*2');
bad_chan = fread(infid,[1,Nbad_chan],'integer*2');
mask = fread(infid,[NEpoch,NChan+1],'integer*2');

