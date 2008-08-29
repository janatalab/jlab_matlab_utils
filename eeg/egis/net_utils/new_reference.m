function [ref_trialdata,mask] = new_reference(trialdata,mask,reference);
%
%ref_trialdata = new_reference(trialdata,mask,reference);
%
% trialdata = trial of data, obviously
% mask = array of zeros and ones with ones indicating good chan
% reference = text or channel_list indicating valid reference type: 
%these types are:
%	'average' = average reference
%	'avgmast' = mathematically averaged mastoids 
%	'perimeter' = perimeter reference (see Shaw thesis for logic)
%	channel_list = list of channels to average together and use as 
%		reference, e.g. [7 86], will reference the data to the average
%		of channels 7 and 86.  
%
if ~(nargin  == 2|nargin == 3)
	error('incorrect number of input arguments')
end

if nargin == 2
	reference = 'average';
end;

if size(trialdata,2)+1 ~= size(mask,2)
	error('mismatch between trialdata and mask');
end;

if strcmp(reference,'average')
	ref_trialdata = average_reference(trialdata,mask);
end;

if strcmp(reference,'vertex')
	ref_trialdata = zeros(size(trialdata,1),size(trialdata,2)+1);
	ref_trialdata(1:size(trialdata,1),1:size(trialdata,2)) = trialdata;
end
if strcmp(reference,'avgmast')
	channel_list = [57 101];
elseif strcmp(reference,'perimeter')
	channel_list = [14 8 1 121 115 101 96 90 83 75 70 64 57 45 39 33 26 22];
end;

if ~(strcmp(reference,'average')|strcmp(reference,'vertex'))
	if size(channel_list,2) > 1	
		cmask = zeros(1,size(mask,2));
		cmask(channel_list) = ones(1,size(channel_list,2));
		cmask = mask.*cmask;
		chan_list = find(cmask);
		ref_trialdata = zeros(size(trialdata,1),size(mask,2));
		ref_trialdata(:,1:size(trialdata,2)) = trialdata;
		if size(chan_list,2) > 1
		ref_sig = (mean(trialdata(:,chan_list)'))';
		ref_trialdata = ref_trialdata - ref_sig*ones(1,size(ref_trialdata,2));
		else
		mask == zeros(1,size(cmask,2));
		end;
	else
		cmask = zeros(1,size(mask,2));
		cmask(channel_list) = ones(1,size(channel_list,2));
		cmask = mask.*cmask;
		chan_list = find(cmask);
		ref_trialdata = zeros(size(trialdata,1),size(mask,2));
		if chan_list == []
			mask = zeros(1,size(mask,2))
		else
			ref_sig = trialdata(:,chan_list);
			ref_trialdata = ref_trialdata - ref_sig*ones(1,size(ref_trialdata,2));
		end;
	end
end;







