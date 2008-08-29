function [Epoch,NEpoch_used,Nbad_chan,bad_chan,avgcsdm] = raw_eeg_csdm(rawfname,add_bad_chan,windows,reference,Epoch,NEpoch,Nbad_chan,bad_chan,mask);
%[Epoch,NEpoch_used,Nbad_chan,bad_chan,avgcsdm] = raw_eeg_csdm(rawfname,
%add_bad_chan,windows,reference,Epoch,NEpoch,Nbad_chan,bad_chan,mask);
%
%this routine calculates the cross spectral density function.  It can be either
% be called serially after dealing with artifact editing, or else itself grab a% preexisting edit code file.   
%
%input arguments
%
% rawfname = NSFragger format raw data file
% add_bad_chan = bad_channels to add (optional), e.g. [8 7 92]
% 
% windows = windows in units of epoch seconds to use in the 
%NSFRagger file starting at the 
%beginning of the file, e.g., [1:6 100:152 170:191] specifies 82 seconds to
%be included in the calculation starting at the beginning of the file, if an
%Epoch length of 1 second is requested, either during editing or within this
%call.  If Epoch =2 this choice of windows implies 164 seconds to be used in
%the average. This parameter defaults to the entire file (optional). 
% 
% reference = reference to use.  The valid types of references are found in
%the m file new_reference.m If desired more should be added there. e.g., 
%'average' to indicate the average reference, defaults to average reference 
%(optional) 'laplacian' can be used as can 'cortical05', 'cortical10', etc. 
%where the number indicates the level of smoothing
%
% the remaining parameters are only used if the edit code information was 
% calculated in the pevious call and is sitting in memory. These parameters 
% are from the '.mask' header.  If you are reading it internally to this       
% m-file do not provide those parameters.  
% (see read_mask.m)
%  	
%output arguments
%
% Epoch = epoch_length in seconds
% NEpoch = actual number of epochs used in spectral calculations
% NBad_chan = actual number of bad_channels
% bad_chan = updated bad channel list
% avgcsdm = cross spectral density function
%
% writtin by R.S. V.1.0
%
% Note: this is being provided to you for your research
% please do not redistribute without my permission
%
if nargout ~= 5
	error('too many or too little output arguments')
end; 
if nargin < 1
	error('filename not specified');
end;
if isempty(rawfname)
	[fid, rawfname] = get_fid('rb');
	fclose(fid)
	fid = fopen(rawfname,'rb','b');
	if nargin <= 4
			[Epoch,NEpoch,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask] = read_mask(rawfname);
	end;		

else
	if (nargin <= 4)
	fid = fopen(rawfname,'rb','b');
	if (fid < 0)
		error('error opening data file')
		fclose('all')
	end;
	[Epoch,NEpoch,Max_microv,Min_Chan,Nbad_chan,bad_chan,mask] = read_mask(rawfname);
	else
	fid = fopen(rawfname,'rb','b');
	if (fid < 0)
		error('error opening data file')
		fclose('all')
	end
	end;
end;	
if ~(nargin ==1 | nargin ==2|nargin == 3| nargin ==4|nargin == 9)
	error('inproper number of input arguments')
end;
if nargin == 1
	add_bad_chan = [];
	windows = [1:size(mask,1)];
	reference = [];
end;
if nargin == 2
	windows = [1:size(mask,1)];
	reference = [];
end;
if nargin == 3
	reference = [];
end;
if isempty(windows)
	windows = [1:size(mask,1)];
end;
bad_chan = [bad_chan add_bad_chan];
Nbad_chan = Nbad_chan + size(add_bad_chan,2);
mask(:,bad_chan) = zeros(size(mask,1),size(bad_chan,2));
mask_win = zeros(size(mask,1),size(mask,2));
if ~isempty(windows)
	mask_win(windows,:) = ones(size(windows,2),size(mask,2));
end;
mask = mask.*mask_win;
[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
avgcsdm = zeros(Samp_Rate*Epoch/2,(NChan+2)*(NChan + 1)/2);
if strcmp(reference,'laplacian')
	[xelec,yelec,zelec] = electrodes(NChan+1);
end;
if strcmp(reference(1:2),'co')
	[xelec,yelec,zelec] = electrodes(NChan+1);
	sources = xyz2tp(xelec,yelec,zelec);
	A = transfer_matrix(500,0.2,80,8,8.2,8.7,9.2,7.8,sources,sources);
end;
for i = 1:NEpoch
if i <= size(mask,1)		
	if sum(mask(i,:)) > 1
		trialdata = fread(fid,[NChan+NEvent,Samp_Rate*Epoch],'integer*2');
		trialdata2 = trialdata(1:NChan,:)';
		trialdata2 = trialdata2*scale;
		if ~isempty(reference)
			if ~(strcmp(reference,'laplacian')|strcmp(reference(1:2),'co'))
			[ref_trialdata,masknew] = new_reference(trialdata2,mask(i,:),reference);
			ref_trialdata = ref_trialdata.*(ones(size(ref_trialdata,1),1)*masknew);
			elseif strcmp(reference,'laplacian')
			good_chan = find(mask(i,:));	
			ref_trialdata = laplacian_trial([trialdata2 zeros(size(trialdata2,1),1)],good_chan,xelec,yelec,zelec);
			masknew = mask(i,:);
			elseif strcmp(reference(1:8),'cortical')					good_chan = find(mask(i,:));
			good_trial = zeros(size(trialdata2,1),size(good_chan,2));
			trialdata3 = average_reference(trialdata2(:,1:NChan),mask(i,:));
			good_trial = trialdata3(:,good_chan);
			good_A = zeros(size(good_chan,2),size(good_chan,2));
			good_A = A(good_chan,good_chan);
			ref_trialdata = zeros(size(trialdata3,1),size(trialdata3,2));
			sigma_m = 1;
			sigma_v = str2num(reference(9:10));
			
			[m,deviations] = bayes_dipole_trial(good_A,good_trial',sigma_v,sigma_m);
			ref_trialdata(:,good_chan) = m';
			masknew = mask(i,:);
			end
		else
			ref_trialdata = average_reference(trialdata2,mask(i,:));
		end;
		
		avgcsdm = avgcsdm + csdm(ref_trialdata);
		mask(i,:) = masknew;
	else
		skip_bytes = (NChan+NEvent)*Samp_Rate*Epoch*2;
		fseek(fid,skip_bytes,'cof');
		
	end;
	
end;
end;
num_good_trials = sum(mask);
icount = 1;
good_cross = zeros(1,size(avgcsdm,2));
for ichan = 1:NChan+1
	for jchan = ichan:NChan+1
		good_cross(icount) = min([num_good_trials(ichan) num_good_trials(jchan)]);
		if good_cross(icount) == 0
			good_cross(icount) = 1;
		end;
		icount = icount+1;
	end;
end;
	
avgcsdm = avgcsdm./(ones(size(avgcsdm,1),1)*good_cross);

NEpoch_used = max(good_cross); 










