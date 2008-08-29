function status = image_meg_power(freqs,obs,dimen,paging,plot_res,power_type,powerfilename,output_name);
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
	fclose(fid);
elseif isempty(powerfname)
	[fid, powerfname, pathname] = get_fid('rb');
	fclose(fid);
end;

pfid = open_file_w_byte_order(powerfname,-2);
[pversion,nfiles,NChan,NMeg,NReference,NEeg,NAnalog,NBad_chan,bad_chan,NEpoch,Epoch,nfreq,obs_labels] = rd_meg_anal_hdr(pfid);

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
maxpower = 0;
figure
iline = 1;
if strcmp(paging,'byobs')
	for i = 1:nfiles
		power = fread(fid,[NChan(i),nfreq(i)],'float');
		if obs_mask(i) > 0 
		power = power';
		if strcmp(power_type,'amplitude')
			power = sqrt(power);
		end;
		bins = round(freqs*Epoch(i) + ones(1,size(freqs,2)));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		imagename = [output_name '.f' int2str(freqs) '.pow'];
		make_meg_image(power(bins,1:NMeg),[57 123 bchan],plot_res,dimen(1),dimen(2),imagename);
		end;
	end;
elseif strcmp(paging,'cross') 
	image_binary = zeros(dimen(1)*plot_res,dimen(2)*plot_res);
	iff = 1;
	for i = 1:nfiles
		power = fread(fid,[NChan(i),nfreq(i)],'float');
		if obs_mask(i) > 0 
		power = power';
		if strcmp(power_type,'amplitude')
			power = sqrt(power);
		end;
		bins = round(freqs*Epoch(i) + ones(1,size(freqs,2)));
		bchan = bad_chan(i,find(bad_chan(i,:)));
		imagename = [output_name '.' deblank(setstr(obs_labels(i,:))) '.f' int2str(freqs)];
		image_binary((iff-1)*plot_res+1:iff*plot_res,:) = make_meg_image(power(bins,1:NMeg),[57 123 bchan],plot_res,1,size(bins,2),[]);
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
	num_written = fwrite(ifid,image_binary','float');
	if num_written ~= size(image_binary,1)*size(image_binary,2)
		error('write-out failed')
	end;
end;	
status = 1;









