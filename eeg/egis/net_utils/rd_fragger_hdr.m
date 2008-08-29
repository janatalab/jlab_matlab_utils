function [header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
%[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
%
% reads in header version 2 of NSFragger output files 
% header_array = complete header except event codes
% EventCodes = event codes 
% Samp_Rate = sampling rate
% NChan = #of channels
% scale = constant to use to convert data to microvolts
% NSamp = sampling rate
% NEvent = number of events
%
% Release 1.0 2/9/97 R.S. 
% 
% Note: This code was developed on my own time, for your lab use.  
%       Please don't redistribute on your own without asking me.  
% 
if (nargin < 1) 
	error('you didnt provide an fid');
end;
if fid < 1
	error('you didnt provide a valid file');
end;
version = fread(fid,1,'integer*4');
if version ~= 2
	error('the jig is up, EGI has updated its software without informing me');
end;

year = fread(fid,1,'integer*2');
month = fread(fid,1,'integer*2');
day = fread(fid,1,'integer*2');
hour = fread(fid,1,'integer*2');
minute = fread(fid,1,'integer*2');
second = fread(fid,1,'integer*2');
millisecond = fread(fid,1,'integer*4');
Samp_Rate = fread(fid,1,'integer*2');
NChan = fread(fid,1,'integer*2');
Gain = fread(fid,1,'integer*2');
Bits = fread(fid,1,'integer*2');
Range = fread(fid,1,'integer*2');
scale = Range/(2^Bits);
NSamp = fread(fid,1,'integer*4');
NEvent = fread(fid,1,'integer*2');
for i = 1:NEvent
	EventCodes(i,1:4) = fread(fid,[1,4],'char*1');
end;
header_array = [version year month day hour minute second millisecond Samp_Rate NChan Gain Bits Range NSamp NEvent];






