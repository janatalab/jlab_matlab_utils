function [ref_trialdata, bad_chan] = grab_ns_data(seconds,reference,Max_microv,rawfname);
%[ref_trialdata, bad_chan] = grab_ns_data(seconds,reference,Max_microv,rawfname);
%
%seconds = seconds from start of file to grab data, e.g., [1:4] or [5] or [97:111]
%reference = reference to use (see new_reference.m) (optional)
%		options: 'average', 'avgmast', 'perimeter', 'vertex',or a channel list
% 		defaults to 'average'
%		can also be used to invoke laplacian: 'laplacian'
%		can also be used to invoke cortical filtering: 'cortical05'
%			or 'cortical10' where the number indicates the
%			smoothing level
%data is alway zeromeaned channel by channel
%Max_microv = maximum microvolt level at channels (optional)
%				defaults to 100
%rawfname = filename (optional)
%
if nargin < 1
	error('seconds not specified');
end;
if (nargin < 2);
	reference = 'average';
end;
if reference == [];
	reference = 'average';
end;
if (nargin < 3);
	Max_microv = 100;
end;
if Max_microv == [];
	Max_microv = 100;
end;
if nargin < 4
	[fid, fname, pathname] = get_fid('rb','*.*');
		rawfname = [pathname fname];
else
	fid = fopen(rawfname,'rb','ieee-be');
end;

[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);

max_seconds = fix(NSamp/Samp_Rate);
if max(seconds) > max_seconds
	fclose('all');
	error('seconds greater than file length');
end;
trialdata = zeros(NChan+NEvent,Samp_Rate*size(seconds,2));	
skip_bytes = Samp_Rate*(seconds(1)-1)*(NChan+NEvent)*2;

fseek(fid,skip_bytes,'cof');
Max_Bins = (Max_microv/scale).^2;
mask = ones(1,NChan+1);
trialdata = fread(fid,[NChan+NEvent,Samp_Rate*size(seconds,2)],'int16');
trialdata2 = zeros(NChan+1,Samp_Rate*size(seconds,2));
trialdata2(1:NChan,:) = trialdata(1:NChan,:).*scale;
trialdata2 = trialdata2';
trialdata(1:NChan,:) = trialdata(1:NChan,:).^2;
max_trial = max(trialdata(1:NChan,:)');
bad_chan = find(max_trial > Max_Bins);
mask(bad_chan) = zeros(1,size(bad_chan,2));
if ~(strcmp(reference,'laplacian')|strcmp(reference(1:2),'co'))
	ref_trialdata = new_reference(trialdata2(:,1:NChan),mask,reference);
elseif strcmp(reference,'laplacian')
	[xelec,yelec,zelec] = electrodes(NChan + 1);
	good_chan = find(mask);
	
	ref_trialdata = laplacian_trial(trialdata2,good_chan,xelec,yelec,zelec);
elseif strcmp(reference(1:8),'cortical')
	[xelec,yelec,zelec] = electrodes(NChan+1);
	sources = xyz2tp(xelec,yelec,zelec);
	A = transfer_matrix(500,0.2,80,8,8.2,8.7,9.2,7.8,sources,sources);
	good_chan = find(mask);
	good_trial = zeros(size(trialdata2,1),size(good_chan,2));
	trialdata3 = average_reference(trialdata2(:,1:NChan),mask);
	good_trial = trialdata3(:,good_chan);
	good_A = zeros(size(good_chan,2),size(good_chan,2));
	good_A = A(good_chan,good_chan);
	ref_trialdata = zeros(size(trialdata3,1),size(trialdata3,2));
	sigma_m = 1;
	sigma_v = str2num(reference(9:10));
	[m,deviations] = bayes_dipole_trial(good_A,good_trial',sigma_v,sigma_m);
	ref_trialdata(:,good_chan) = m';
end;
ref_trialdata = zeromean(ref_trialdata);


