function status = lap_ns_data(rawfname,Max_microv,Min_Chan);
%status = lap_ns_data(rawfname);
%rawfname = filename (optional)
%Max_microv = maximum micovolt value at a channel (optional)
%				defaults to 100
%Min_Chan = minimum number of channels to use trial (optional)
%				defaults to 100
%takes the surface Laplacian of a NSFragger file creating new
%NSFragger file.  The new file has NChan + 1 channels.  Event 
%codes are transmitted to the new file.  
if nargin < 1
		[fid, fname, pathname] = get_fid('rb','*.*');
		rawfname = [fname];
else
	    fid = fopen(rawfname,'rb');
end;

if rawfname == []
	[fid, fname, pathname] = get_fid('rb','*.*');
end
if nargin < 2
	Max_microv = 100;
	Min_Chan = 100;
end;

[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
[Epoch, NEpoch,Max_microv,Min_Chan,NBad_chan,bad_chan,mask] = ns_artifact_edit(fid,Samp_Rate,NChan, NSamp, NEvent,scale,1,Max_microv,Min_Chan);
frewind(fid);
[header_array, EventCodes,Samp_Rate, NChan, scale, NSamp, NEvent] = rd_fragger_hdr(fid);
header_array(13) = header_array(13)/100;
header_array(10) = NChan + 1;
header_array(14) = Samp_Rate*NEpoch;
[xelec,yelec,zelec] = electrodes(NChan+1);
lap_filename = [rawfname '_lap'];
lapfid = fopen(lap_filename,'wb');

wt_fragger_hdr(lapfid,header_array,EventCodes);
for i = 1:NEpoch
	if sum(mask(i,:)) > 1
		trialdata = zeros(NChan+NEvent,Samp_Rate);
		trialdata = fread(fid,[NChan+NEvent,Samp_Rate],'int16');
		trialdata = trialdata';
		trialdata2 = zeros(Samp_Rate,NChan+1);
		trialdata2(:,1:NChan) = trialdata(:,1:NChan);
		good_chan = find(mask(i,:));
		lap_trialdata2 = laplacian_trial(trialdata2,good_chan,xelec,yelec,zelec);
	    plot(lap_trialdata2(:,73));
		lap_trialdata2 = 100*lap_trialdata2';
		lap_trialdata = zeros(NChan+1+NEvent,Samp_Rate);
		lap_trialdata(NChan+2:NChan+1+NEvent,:) = trialdata(NChan+1:NChan+NEvent,:);
		lap_trialdata(1:NChan+1,:) = lap_trialdata2;
		figure
		plot(lap_trialdata(73,:));
	else
		lap_trialdata = zeros(NChan+1+NEvent,Samp_Rate);
	end
	fwrite(lapfid,lap_trialdata,'int16');
end;
fclose(fid);
fclose(lapfid);
	
