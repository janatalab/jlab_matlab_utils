function status = run_eeg_csdm(rawfname,add_bad_chan,windows,reference,max_freq,Epoch,Max_microv,Min_Chan,LoPass,HiPass);
%status = run_eeg_csdm(rawfname,add_bad_chan,windows,reference,max_freq,
%			Epoch,Max_microv,Min_Chan,LoPass,HiPass);
%
%makes cross spectral density function files (version = -1 file)
%from a NSFragger raw2  file. This is the main routine for spectral analysis
%
%rawfname = cross spectral density filename. Can submit a blank [] for gui
%add_bad_chan = additional bad channels over artifact editing (optional)
%windows = epochs to use in the analysis.  By epochs I mean that the units 
%	of this parameter are epochs which are either specified in the edit 
%	code file or as the parameter Epoch, e.g. [2:4 7], means use epochs
%	2,3,4 and 7 in the average.  If Epoch = 1 this is 4 seconds of EEG;
%	if Epoch = 2 this is 8 seconds of EEG but still only 4 epochs in the 
%	average (optional - default is whole file)
%reference = reference to use legal values are: 'average', 'avgmast', 
%	'perimeter','laplacian','vertex'. Or else a channel list can be
%	submitted, e.g. [5 6 7], or 'cortical05', 'cortical10' for cortical
%	filtering routines. 
%	optional - default is 'average' 
%max_freq = highest frequency to save data (must be less than SampRate/2)
%			default is 50 Hz
%If you have previously edited the file stop here. if not add the following 
%params to edit 
%Epoch = epoch length to use 
%Max_microv = maximumvoltage in microvolts
%Min_Chan = minimum number of channels
%LoPass = optional lowpass filter settings (optional)
%HiPass = optional hipass settings. (optional)
if nargin < 1
	error('filename not specified');
end;
if isempty(rawfname)
	[fid, rawfname] = get_fid('rb');
	if (fid < 0)
		error('error opening data file')
		fclose('all')
	end;
	fclose(fid);
	fid = fopen(rawfname,'rb','b');
else
	fid = fopen(rawfname,'rb','b');
	if (fid < 0)
		error('error opening data file')
		fclose('all')
	end;
end;

if ~(nargin ==1 | nargin ==2|nargin == 3| nargin ==4|nargin == 5|nargin == 9|nargin == 10|nargin == 8)
	error('inproper number of input arguments')
end;
if nargin < 8
	[Epoch,NEpoch,Max_microv,Min_Chan,Nbad_chan,bad_chan,mmask] = read_mask(rawfname);
end;
if nargin == 1
	add_bad_chan = [];
	windows = [1:size(mmask,1)];
	reference = ['average'];
	max_freq = 50;
end;
if nargin == 2
	windows = [1:size(mmask,1)];
	reference = ['average'];
	max_freq = 50;
end;
if nargin == 3
	reference = ['average'];
	max_freq = 50;
end;
if nargin == 4
	max_freq = 50;
end;

if (isempty(windows) & nargin < 8)
	windows = [1:size(mmask,1)];
end;
if isempty(reference)
	reference = 'average';
end;
if isempty(max_freq)
	max_freq = 50;
end

if nargin == 10
	[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvenp] = rd_fragger_hdr(fid);
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mmask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass,Hi_Pass); 
end
if nargin == 9
	[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mmask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan,Lo_Pass);
end;
if nargin == 8
	[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
	[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mmask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,Epoch,Max_microv,Min_Chan);
end;
fclose(fid);
if isempty(windows)
	windows = [1:size(mmask,2)];
end;
for i = 1:size(windows,1)
	if size(windows,1) > 1	
		if i < 10
			text_lab = ['.' reference(1,1:4) '.w.0'];
		else
			text_lab = ['.' reference(1,1:4) '.w.'];
		end;
		csdmfname = [rawfname(1:size(rawfname,2)-8) '.csdm.' text_lab int2str(i)];
	else
		csdmfname = [rawfname(1:size(rawfname,2)-8) '.' reference(1,1:4) '.csdm'];
	end;
	csdmfid = fopen(csdmfname,'wb');
	if nargin >= 8
		[Epoch,NEpoch_used,Nbad_chan,bad_chan,avgcsdm] = raw_eeg_csdm(rawfname,add_bad_chan,windows(i,:),reference,Epoch,NEpoch,Nbad_chan,bad_chan,mmask);
	else
		[Epoch,NEpoch_used,Nbad_chan,bad_chan,avgcsdm] = raw_eeg_csdm(rawfname,add_bad_chan,windows(i,:),reference);
	end;
	version = -1;
	fwrite(csdmfid,version,'integer*2');
	fwrite(csdmfid,Epoch,'integer*2');
	fwrite(csdmfid,size(windows,2)*Epoch,'integer*2')
	fwrite(csdmfid,NEpoch_used,'integer*2');	
	fwrite(csdmfid,Nbad_chan,'integer*2');
	fwrite(csdmfid,bad_chan,'integer*2');
	if strcmp(reference,'average')
		ref_flag = 1;
	elseif strcmp(reference,'avgmast')
		ref_flag = 2;
	elseif strcmp(reference,'perimeter');
		ref_flag = 3;
	elseif strcmp(reference,'vertex');
		ref_flag = 4;
	elseif strcmp(reference,'laplacian');
		ref_flag = 5;
	elseif strcmp(reference(1:8),'cortical');
		ref_flag = 7;
	else
		ref_flag = 6;
	end;
	fwrite(csdmfid,ref_flag,'integer*2');
	if ref_flag == 6
		fwrite(csdmfid, size(reference,2),'integer*2');
		fwrite(csdmfid, reference,'integer*2');
	end;
	if ref_flag == 7
		sigma_v = str2num(reference(9:10));
		fwrite(csdmfid,sigma_v,'integer*2');
	end;
	max_freq = max_freq*Epoch+1;
	fwrite(csdmfid,max_freq,'integer*2');
	fwrite(csdmfid,size(avgcsdm,2),'integer*2');
	fwrite(csdmfid,real(avgcsdm(1:max_freq,:))','real*4');
	fwrite(csdmfid,imag(avgcsdm(1:max_freq,:))','real*4');
	fclose(csdmfid);
end;
status = 1;
fclose('all');
	















