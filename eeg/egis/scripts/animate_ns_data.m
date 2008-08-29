function status = animate_ns_data(seconds,reference,rawfname,plot_res,Max_microv);
%status = animate_ns_data(seconds,reference,rawfname,plot_res,Max_microv);
%
%seconds = seconds to animate in NSFragger file
%reference = reference to use in animating the data (optional)
%	'average', 'avgmast', 'perimeter', 'vertex', or a channel list
%	'laplacian' = surface Laplacian estimate
%	defaalts to 'average'
%rawfname = filename (optional)
%Max_microv = maximum microvolts (optional)
%plot_res = size of output plots
if nargin < 1
	error('seconds not specified')
end

if nargin < 3
	if nargin < 2
	reference = 'average';
	[fid, rawfname] = get_fid('rb');
	fclose(fid);
	Max_microv = 100;
	plot_res = 40;
	else
	[fid, rawfname] = get_fid('rb')
	fclose(fid);
	Max_microv = 100;
	plot_res = 40;
	end;
end
if nargin < 4
	plot_res = 40;
	Max_microv = 100;
end
if nargin < 5
	plot_res = 40;
end;
if rawfname == []
	[fid, rawfname] = get_fid('rb')
	fclose(fid);
end;
if reference == []
	reference = 'average';
end;
if Max_microv == []
	Max_microv = 100;
end;
if plot_res == []
	plot_res = 40;
end;

[ref_trialdata, bad_chan] = grab_ns_data(seconds,reference,Max_microv,rawfname);

[xelec,yelec,zelec] = electrodes(size(ref_trialdata,2));
good_chan_mask = ones(1,size(ref_trialdata,2));
ruu = sqrt(xelec(1).^2+yelec(1).^2+zelec(1).^2);
[xs,ys,zs] = polgrid(plot_res,ruu);
if bad_chan ~= []	
good_chan_mask(1,bad_chan) = zeros(1,size(bad_chan,2));
end;
good_chan = find(good_chan_mask(1,:));
x = zeros(1,size(good_chan,2));
y = zeros(1,size(good_chan,2));
z = zeros(1,size(good_chan,2));
x = xelec(good_chan);
y = yelec(good_chan);
z = zelec(good_chan);
welec = 1;
icount = 1;
[k,kinv,a,ainv,e] = k_and_e(welec,x,y,z);
for i = 1:size(ref_trialdata,1)
       	w = ref_trialdata(i,:);
	    v = zeros(1,size(good_chan,2));
		v = w(good_chan);
		[p,q,error_check]= mateqs(welec,x,y,z,v,k,kinv,a,ainv,e);
		image_3d= interp_3d(welec,x,y,z,xs,ys,zs,p,q);	
		if ~strcmp(reference,'laplacian')
			imagemat= mask(xs,ys,zs,image_3d,120); 
		else
			imagemat = mask(xs,ys,zs,image_3d,108);
		end;
		if i > 9
	 	outfname = [rawfname '.img_s' int2str(i) ];
		else
		outfname = [rawfname '.img_s0' int2str(i) ];
		end;
  	  	outfid = fopen(outfname,'wb');
	   	num_written = fwrite(outfid,imagemat','float');
		fclose(outfid);
end;

