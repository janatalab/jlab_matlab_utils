function status = make_eeg_coherence(csdmfname,powerfname,obs_labels,add_bad);
%status = make_eeg_coherence(csdmfname,coherencefname,obs_labels,add_bad);
%generate multisubject/multistate coherence file from multiple csdm
%files
%csdmfname = root (everything short of the numbers at the end) of the csdm
%	file or the entire csdmfilename or csdmfilenames
%	if you plan to use the indices option then what you
%	need to do is have the root without the indice which must be 
%	then appended. 
%	you can only do a single file if you use gui by submitting a blank
%	or no parameters
%	you can submit a list of csdm file names as rows in csdmfname.  In 
%	this case all the csdm filenames must have the same length.
%	default is a gui for a single csdm file  (optional)
%coherence_fname = name of coherence file.
%obs_labels = labels for obbservations (optional),default = csdmfname 
%add_bad = even more bad_channels
%example: to process the csdm files 
%	010_12.EEG.RAW.csdm.aver.w01
%	010_12.EEG.RAW.csdm.aver.w02
%	010_12.EEG.RAW.csdm.aver.w10
%
%or you could list all the csdm files
%csdmfname = [010_12.EEG.RAW.csdm.aver.w01;010_12.EEG.RAW.csdm.aver.w02;
%             010_12.EEG.RAW.csdm.aver.w10];
%make_eeg_coherence(csdmfname)
%or
%make_eeg_coherence(csdmfname,'010_12.coh',[Obs_01;Obs_02;Obs_10]);
%
if nargin < 1
	[fid, csdmfname] = get_fid('rb');
	if fid < 0
		error('invalid filename')
	end;
	fclose(fid);
end;

if nargin < 2
	[pfid,powerfname] = put_fid('rb');
	obs_labels = csdmfname;
	fclose(pfid);
end;

if nargin < 3
	obs_labels = csdmfname;
end;
if isempty(powerfname)
	[pfid,powerfname] = put_fid('rb');
	fclose(pfid);
end;
version = -2;
pfid = fopen(powerfname,'wb');
fwrite(pfid,version,'int16');
fwrite(pfid,size(csdmfname,1),'int16');

fwrite(pfid,size(obs_labels,2),'int16');
for i = 1:size(obs_labels,1)
	fwrite(pfid,obs_labels(i,:),'char');
end;
nfiles = size(csdmfname,1);


PEpoch = zeros(1,nfiles);
PWindow_Length = zeros(1,nfiles);
PNEpoch = zeros(1,nfiles);
PNbad_chan = zeros(1,nfiles);
Pbad_chan = zeros(nfiles,129);
Pref_flag = zeros(1,nfiles);
PNchan = zeros(1,nfiles);
PNfreq = zeros(1,nfiles); 
Preference = zeros(nfiles,129);

	for i = 1:size(csdmfname,1)
		fid = fopen(csdmfname(i,:),'rb');
		frewind(fid);
		version = fread(fid,1,'int16');
		if version ~= -1
			error(['file: ' csdmfname(i,:) 'is not a csdm file']);
			fclose(fid)
		end;
		[PEpoch(i),PWindow_Length(i),PNEpoch(i),PNbad_chan(i),bad_chan,Pref_flag(i),reference,PNchan(i),PNfreq(i)] = rd_csdm_hdr(fid);
		
		if Pref_flag(i) ~= 6
		Preference(i,1:size(reference,2)) = reference;
		else
		Preference(i,1:size(reference,2)) = int2str(reference);
		end;
		if strcmp(reference,'laplacian')
			bad_chan = [bad_chan 22 26 33 39 45 57 64 70 75 83 90 96 101 115 121 1 8 14];
		end
		if nargin == 4
		bad_chan = [bad_chan add_bad(i,:)];
		end;
		Pbad_chan(i,1:size(bad_chan,2)) = bad_chan;	
		PNBad_chan(i) = size(bad_chan,2);
		fclose(fid);
	end;

fwrite(pfid,PEpoch,'int16');
fwrite(pfid,PWindow_Length,'int16')
fwrite(pfid,PNEpoch,'int16');	
fwrite(pfid,PNbad_chan,'int16');
fwrite(pfid,Pbad_chan,'int16');
fwrite(pfid,Pref_flag,'int16');
fwrite(pfid,Preference,'char');
PNchan = (PNchan.*(PNchan+ones(1,nfiles)))/2;		
fwrite(pfid,PNchan,'int16');
fwrite(pfid,PNfreq,'int16');

	for i = 1:size(csdmfname,1)
		fid = open_file_w_byte_order(csdmfname(i,:),-1);
		version = fread(fid,1,'int16');
		if version ~= -1
			error(['file: ' csdmfname(i,:) 'is not a csdm file']);
			fclose(fid)
		end;
		[Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,reference,ravgcsdm,iavgcsdm] = rd_csdm(fid);
		avg_coherence = eeg_coherence(ravgcsdm,iavgcsdm,bad_chan);
		fwrite(pfid,avg_coherence','float');
		fclose(fid);
	end;

fclose(pfid)
status = 1;
	
	



