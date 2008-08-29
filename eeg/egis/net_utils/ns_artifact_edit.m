 function [Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass,Hi_Pass); 
% [Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,
%       Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass,Hi_Pass);
% 
% artifact edits NSFragger file generating mask
%
% fid = legal fid
% Samp_Rate = sampling rate
% NChan = # of channels
% NEvent = # of event tracks
% scale = scale factor between bins and microvolts
% Epoch = analysis epoch length
% Max_microv = ceiling in microvolts at a single channel
% Min_Chan = minimum number of channels to have a legal trial
% Lo_Pass (optional) = low-pass filter frequency  of trial before atifact editing
% Hi_Pass (optional) = hi-pass filter frequency of trial before artifact editing 
%
% Version 1.0 2/8/97 R.S.
%
% Note: The code is provided to you for lab use only
%	please do not retransmit without permission
%

if nargin < 9
	error('not enough parameters')
end;
if fid < 1
	error('invalid fid');
end;
Max_Bins = (Max_microv/scale).^2;

NEpoch = fix(NSamp/(Epoch*Samp_Rate));

mask = ones(NEpoch, NChan+1);
 
if nargin >= 10

[blow,alow] = butter(20,[Lo_Pass/(Samp_Rate/2)]);

end;
if nargin == 11

[bband,aband] = butter(3,[Hi_Pass/(Samp_Rate/2) Lo_Pass/(Samp_Rate/2)]);

end;

for i = 1:NEpoch
	 
	trialdata = fread(fid,[NChan+NEvent,Samp_Rate*Epoch],'int16');
	trialdata = trialdata'; 
	if nargin == 10
	
		for ch =1:NChan
			trialdata(:,ch) = filtfilt(blow,alow,trialdata(:,ch));
		end;
	end;
	if nargin == 11
		for ch =1:NChan
			trialdata(:,ch) = filtfilt(blow,alow,trialdata(:,ch));
			trialdata(:,ch) = filtfilt(bband,aband,trialdata(:,ch));
		end;
	end
	trialdata = trialdata.^2;
	max_trial = max(trialdata);
	bad_chan = find(max_trial > Max_Bins);
	mask(i,bad_chan) = zeros(1,size(bad_chan,2));
	if (sum(mask(i,:)) < Min_Chan)
		mask(i,:) = zeros(1,NChan+1);
	end;
end;

grand_mask = sum(mask);
max_mask = max(grand_mask);
% bad_channels are those which are not good for 75% of trials

bad_chan = find(grand_mask < 0.75*max_mask);
if bad_chan == []
	NBad_chan = 0;
else
	NBad_chan = size(bad_chan,2);
end






