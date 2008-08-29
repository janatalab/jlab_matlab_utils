function status = run_ns_edit(Epoch,Max_microv,Min_Chan,rawfname,Lo_Pass,Hi_Pass);
% status = run_ns_edit(Epoch,Max_microv,Min_Chan,rawfname,Lo_Pass,Hi_Pass);
%
% Epoch = epoch length in seconds for artifact editing and spectral analysis
% Max_microv = maximum value in microvolts at a channel
% Min_Chan = minimum number of channelsin a good trial
% rawfname (optional) = NS Fragger filename
% Lo_Pass (optional) = Lo_Pass filtering frequency
% Hi_Pass (optional) = Hi_Pass filtering frequency
%
% Note: this code was developed for use in this lab
%       please do not redistribute without my permission
%
% Version 1.0 R.S. 2/10/97

if nargin < 3
	error('not enough input arguments');
end;
if nargin == 3;
	[infid, rawfname] = get_fid('rb');
	fclose(infid);
end
if rawfname == [];
	[infid, rawfname] = get_fid('rb');
	fclose(infid);
	infid = fopen(rawfname,'rb','b');
else
	infid = fopen(rawfname,'rb','b');
end;
if infid < 1
	error('invalid filename')
end
outfname = [rawfname '.mask'];
outfid = fopen(outfname,'wb');
fid = infid;
[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(infid);
if nargin == 6
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass,Hi_Pass);
elseif nargin == 5
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass);
else
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan);
end;


status = write_mask(outfid,Epoch,NEpoch,NChan,Max_microv,Min_Chan,NBad_chan,bad_chan,mask);

disp('BAD CHANNELS ARE')
bad_chan

disp('NUMBER OF AVAILABLE TRIALS')
sum(mask(:,NChan+1))










