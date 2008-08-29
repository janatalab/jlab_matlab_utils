function export_matrix = image_eeg_power(freqs,dimen,obs,paging,plot_res,power_type,powerfname,output_name);
%status = = image_eeg_power(freqs,dimen,obs,paging,plot_res,power_type,powerfname,output_name);
%
%freqs = frequencies to plot, e.g. [3:7 10] is 3 4 5 6 and 7 and 10 
% [3:0.5:5] is 3 , 3.5 4,4.5,5, default is 5:1:20 
%dimen = dimensions of image e.g., [2 3] default is [sqrt(size(freqs,2)] sqrt(size(freqs,2) 
%obs = observations (default is all)
%paging = 'byobs' or 'cross' default is 'byobs', cross means frequencies
%	are columns and observations are rows. 
%	if 'cross' is selected the output is 
%	image.o (observation nums) .f (freqs).img
%	if 'byobs' is selected the output is 
%	obs_label.f (freqs) .img
%plot_res = plot resolution e.g., 80 default is 50
%power_type = 'power' or 'amplitude' (default)
%powerfname = powerfilename (can be gui by skipping or passing a blank)
%output_name = root name of output file, defaults to gui
if nargin < 7
	[fid, powerfname, pathname] = get_fid('rb');
	output_name = powerfname;
elseif isempty(powerfname)
	[fid, powerfname, pathname] = get_fid('rb');
else
	fid = fopen(powerfname,'rb');
	pathname = [];
end;
version = fread(fid,1,'int16');
if version ~= -3
	error('this is not a power file');
end;
 [nfiles,obs_labels,Epoch,Window_Length,NEpoch,Nbad_chan,bad_chan,ref_flag,reference,NChan,NFreq] = rd_anal_hdr(fid);
if nargin == 0
	freqs = [5:1:20];
	obs = [1:nfiles];
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
	paging = 'byobs';
	plot_res = 50;
	power_type = 'amplitude';
end;
if nargin == 1
	obs = [1:nfiles];
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
	paging = 'byobs';
	plot_res = 50;
	power_type = 'amplitude';
end;		
if nargin == 2
	obs = [1:nfiles];
	paging = 'byobs';
	plot_res = 50;
	power_type = 'amplitude';
end;	
if nargin == 3
	paging = 'byobs';
	plot_res = 50;
	power_type = 'amplitude';
end;
if nargin == 4
	plot_res = 50;
	power_type = 'amplitude';
end;
if nargin == 5
	power_type = 'amplitude';
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
if isempty(power_type)
	power_type = 'amplitude';
end;
if nargin < 8
	[outfid,output_name] = put_fid('wb')
	fclose(outfid);
end
if isempty(output_name)
	output_name = powerfname;
end;
obs_mask = zeros(1,nfiles);
obs_mask(obs) = ones(1,size(obs,2));
if strcmp(paging,'byobs')
	for i = 1:nfiles
		power = fread(fid,[NChan(i),NFreq(i)],'float');
		if obs_mask(i) > 0 
		power = power';
		if strcmp(power_type,'amplitude')
			power = sqrt(power);
		end;
		bins = freqs*Epoch(i) + ones(1,size(freqs,2));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		imagename = [output_name '.f' int2str(freqs) '.pow'];
		export_matrix = make_a_image(power(bins,:),bchan,plot_res,dimen(1),dimen(2),imagename);
		end;
	end;
elseif strcmp(paging,'list')
	if size(freqs,2) ~= size(obs,2)
		error('lists dont match');
	end;
	if isempty(dimen)
	dimen = [sqrt(size(freqs,2)) sqrt(size(freqs,2))];
	end;
	icount = 1;
	bchan = [];
	for i = 1:nfiles
		power = fread(fid,[NChan(i),NFreq(i)],'float');
		if obs_mask(i) > 0 
		power = power';
		if strcmp(power_type,'amplitude')
			power = sqrt(power);
		end;
		bin = freqs(icount)*Epoch(i) + 1;
		power_plot(icount,:) = power(bin,:);
		bchan = [bchan bad_chan(i,find(bad_chan(i,:)))];
		icount = icount + 1;	
		end;
	end;
		imagename = [output_name];
		export_matrix = make_a_image(power_plot,bchan,plot_res,dimen(1),dimen(2),imagename);
	
elseif strcmp(paging,'cross') 
	image_binary = zeros(dimen(1)*plot_res,dimen(2)*plot_res);
	iff = 1;
	for i = 1:nfiles
		power = fread(fid,[NChan(i),NFreq(i)],'float');
		if obs_mask(i) > 0 
		power = power';
		if strcmp(power_type,'amplitude')
			power = sqrt(power);
		end;
		bins = freqs*Epoch(i) + ones(1,size(freqs,2));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		imagename = [output_name '.' deblank(setstr(obs_labels(i,:))) '.f' int2str(freqs)];
		image_binary((iff-1)*plot_res+1:iff*plot_res,:) = make_a_image(power(bins,:),bchan,plot_res,1,size(bins,2),[]);
		iff = iff + 1;
		end;
	end;
	imagename = [output_name '.o_' int2str(obs) '.f_' int2str(freqs) '.pow.img'];
	ifid = fopen(imagename,'wb');
	if ifid < 0
		fclose('all');
		error('file open for image failed');
	end;
	num_written = fwrite(ifid, size(image_binary), 'float');
	if num_written ~= 2
			error('Failed to write header info');
	end 
	fwrite(ifid,plot_res,'float');
	export_matrix = image_binary;
	num_written = fwrite(ifid,image_binary','float');
	if num_written ~= size(image_binary,1)*size(image_binary,2)
		error('write-out failed')
	end;
end;	
status = 1;








