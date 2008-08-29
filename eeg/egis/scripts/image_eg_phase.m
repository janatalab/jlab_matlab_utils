function status = image_eeg_coherence(freqs,dimen,obs,paging,plot_res,ref_chan,coherencefname);
%status = = image_eeg_coherence(freqs,dimen,obs,paging,plot_res,ref_chan,coherencefname);
%
%freqs = frequencies to plot, e.g. [3:7 10] is 3 4 5 6 and 7 and 10 
% [3:0.5:5] is 3 , 3.5 4,4.5,5, default is 5:1:20 
%dimen = dimensions of image e.g., [2 3] default is [sqrt(size(freqs,2)] sqrt(size(freqs,2) 
%obs = observations (default is all)
%paging = 'byobs' or 'cross' default is 'byobs', cross means frequencies
%	are columns and observations are rows. 
%	if 'cross' is selected the output is 
%	image.o (observation nums) .f (freqs).coh.img
%	if 'byobs' is selected the output is 
%	obs_label.f (freqs) .coh.img
%plot_res = plot resolution e.g., 80 default is 50
%ref_chan = channel to plot coherence w/respect to (default is 73)
%coherencefname = coherencefilename (can be gui by skipping)
%
if nargin < 8
	[fid, fname, pathname] = get_fid('rb');
else
	fid = fopen(coherencefname,'rb');
	pathname = [];
end;
version = fread(fid,1,'int16');
if version ~= -2
	error('this is not a coherence file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [5:1:20];
	obs = [1:nfiles];
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
	paging = 'byobs';
	plot_res = 50;
	ref_chan = 73;
end;
if nargin == 1
	obs = [1:nfiles];
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
	paging = 'byobs';
	plot_res = 50;
	ref_chan = 73;
end;		
if nargin == 2
	obs = [1:nfiles];
	paging = 'byobs';
	plot_res = 50;
	ref_chan = 73;
end;	
if nargin == 3
	paging = 'byobs';
	plot_res = 50;
	ref_chan = 73;
end;
if nargin == 4
	plot_res = 50;
	ref_chan = 73;
end;
if nargin == 5
	ref_chan = 73;
end;
if isempty(paging)
	paging = 'byobs';
end
if isempty(plot_res)
	plot_res = 50;
end;
if isempty(dimen)
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
end

if isempty(obs)
	obs = [1:nfiles];
end;
if isempty(ref_chan)
	ref_chan = 73;
end;

obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
ch_pair_indices;
if strcmp(paging,'byobs')
	for i = 1:nfiles
		if obs_mask(i) > 0 
		avgcoherence = fread(fid,[NChan(i),NFreq(i)],'float');
		avgcoherence = avgcoherence';
		bins = freqs*Epoch(i) + ones(1,size(freqs,2));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		imagename = [pathname deblank(setstr(obs_labels(i,:))) '.f' int2str(freqs) '.coh'];
		pchan = chpair2nchan(NChan(i));
		cohplt = zeros(bins,pchan);
		for i = 1:pchan
			cohplt(bins,i) = avgcoherence(bins,chpair(ref_chan,i));
		end;
		make_a_image(cohplt(bins,:),bchan,plot_res,dimen(1),dimen(2),imagename);
		else
		skip_bytes = NChan(i)*NFreq(i)*4;
		fseek(fid,skip_bytes,'cof');
		end;
	end;
elseif strcmp(paging,'cross') 
	image_binary = zeros(dimen(1)*plot_res,dimen(2)*plot_res);
	iff = 1;
	for i = 1:nfiles
		if obs_mask(i) > 0 
		avgcoherence = fread(fid,[NChan(i),NFreq(i)],'float');
		avgcoherence = avgcoherence';
		bins = freqs*Epoch(i) + ones(1,size(freqs,2));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		pchan = chpair2nchan(NChan(i));
		for i = 1:pchan
			cohplt(bins,i) = avgcoherence(bins,chpair(ref_chan,i));
		end;
		image_binary((iff-1)*plot_res+1:iff*plot_res,:) = make_a_image(cohplt(bins,:),bchan,plot_res,1,size(bins,2),[]);
		iff = iff + 1;
		else
		skip_bytes = NChan(i)*NFreq(i)*4;
		fseek(fid,skip_bytes,'cof');
		end;
	end;
	imagename = [pathname 'image.o' int2str(obs) '.f' int2str(freqs) '.coh.img'];
	ifid = fopen(imagename,'wb');
	if ifid < 0
		fclose('all');
		error('file open for image failed');
	end;
	num_written = fwrite(ifid, size(image_binary), 'float');
	if num_written ~= 2
			error('Failed to write header info');
	end 

	num_written = fwrite(ifid,image_binary','float');
	if num_written ~= size(image_binary,1)*size(image_binary,2)
		error('write-out failed')
	end;
end;	
status = 1;








