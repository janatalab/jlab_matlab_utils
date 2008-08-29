function status = write_mask(fid,Epoch,NEpoch,NChan,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask);
%
% status = write_mask(fid,Epoch,NEpoch,NChan,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask);
%
%writes artifact editing information to disk for use in other analysis
%
%
if nargin < 8 
	error('insufficient input arguments')
end;

if fid < 1
	error('invalid file id')
end;

fwrite(fid,Epoch,'int16');
fwrite(fid,NEpoch,'int16');
fwrite(fid,NChan,'int16');
fwrite(fid,Max_microv,'int16');
fwrite(fid,Min_Chan,'int16');
fwrite(fid,Nbad_chan,'int16');
fwrite(fid,bad_chan,'int16');
fwrite(fid,mask,'int16');
status = 1;
