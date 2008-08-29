function status = wt_fragger_hdr(fid,header_array,EventCodes);
%status = wt_fragger_hdr(fid,header_array,EventCodes);
%
%write a header to a new NSFragger file
%
%fid = open new file
%header_array = header
%EventCodes = event codes
% 
% Release 1.0 2/9/97 R.S. 
% 
% Note: This code was developed on my own time, for your lab use.  
%       Please don't redistribute on your own without asking me.  
% 
if nargin < 3
 	error('insufficient input arguments');
end
if fid < 1
	error('you have provided an invalid file');
end;
if size(header_array,2) ~= 15
	error('size of header array is incorrect');
end;
fwrite(fid,header_array(1),'int32');
fwrite(fid,header_array(2),'int16');
fwrite(fid,header_array(3),'int16');
fwrite(fid,header_array(4),'int16');
fwrite(fid,header_array(5),'int16');
fwrite(fid,header_array(6),'int16');
fwrite(fid,header_array(7),'int32');
fwrite(fid,header_array(8),'int16');
fwrite(fid,header_array(9),'int16');
fwrite(fid,header_array(10),'int16');
fwrite(fid,header_array(11),'int16');
fwrite(fid,header_array(12),'int16');
fwrite(fid,header_array(13),'int16');
fwrite(fid,header_array(14),'int32');
fwrite(fid,header_array(15),'int16');
for i = 1:header_array(15)
	fread(fid,EventCodes(i),'char');
end;
status =1;

