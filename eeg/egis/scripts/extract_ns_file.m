function status= extract_ns_file(seconds,rawfname);
%status= extract_ns_file(seconds,rawfname);
%seconds = seconds from start of file to grab data, e.g., [1:4] or
%  [5] or [97:111]
%rawfname = filename (optional)
%
if nargin < 1
	error('seconds not specified');
end;

if nargin < 2
	[fid, fname, pathname] = get_fid('rb','*.*');
		rawfname = [fname];
else
	fid = fopen(rawfname,'rb');
end;

outfname = [rawfname '_ext_' int2str(min(seconds)) '_' int2str(max(seconds))]
outfid = fopen(outfname,'wb');
[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
header_array(14) = size(seconds,2)*Samp_Rate;
wt_fragger_hdr(outfid,header_array,EventCodes);
max_seconds = fix(NSamp/Samp_Rate);
if max(seconds) > max_seconds
	fclose('all');
	error('seconds greater than file length');
end;
trialdata = zeros(NChan+NEvent,Samp_Rate);
	
skip_bytes = Samp_Rate*(seconds(1)-1)*(NChan+NEvent)*2;

fseek(fid,skip_bytes,'cof');
for i = 1:size(seconds,2)
trialdata = fread(fid,[NChan+NEvent,Samp_Rate],'int16');
	fwrite(outfid,trialdata,'int16');
end
fclose('all');
status = 1;
