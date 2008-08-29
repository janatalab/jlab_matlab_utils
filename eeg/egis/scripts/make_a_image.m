function [image_binary, max_mask] = make_a_image(data,bad_chan,plot_res,image_rows,image_columns,out_root);
%status = make_a_image(data,bad_chan,plot_res,image_rows,image_columns
%,out_root)
%data = arbitrary number of rows Nchan columns. 
%bad_chan = bad_Chans list e.g [6 7 9]; (optional)
%plot_Res = image size e.g 100 (optional) defaults to 100
%e.g.
%make_a_image(cor_by_obs([1 4],:),[37 68],100)
%other optional parameters
%image_rows (default is 1)
%image_columns (default is size(data,1))
%out_root (default is gui) generates file out_root.img
%
%note: thisis the lowest level imaging software, there is probably a higher
%level script you should be using
if ~ (nargin ==1|nargin ==2|nargin == 3|nargin == 5|nargin == 6)
	error('inappropriate number of input arguments');
end;
 
if nargin < 2
	bad_chan = [];
	plot_res = 100;
end
if nargin < 3
	plot_res= 100;
end;
if nargin <  4
image_binary = zeros(plot_res,plot_res*size(data,1));
image_rows = 1;
image_columns = size(data,1);
else
image_binary = zeros(plot_res*image_rows,plot_res*image_columns);
end;

[xelec,yelec,zelec] = electrodes(size(data,2));
good_chan_mask = ones(1,size(data,2));

ruu = sqrt(xelec(1).^2+yelec(1).^2+zelec(1).^2);
[xs,ys,zs] = polgrid(plot_res,ruu);
if ~isempty(bad_chan)	
good_chan_mask(1,bad_chan) = zeros(1,size(bad_chan,2));
end;
good_chan = find(good_chan_mask(1,:));
x = zeros(1,size(good_chan,2));
y = zeros(1,size(good_chan,2));
z = zeros(1,size(good_chan,2));
x = xelec(good_chan);
y = yelec(good_chan);
z = zelec(good_chan);
elec = xyz2tp(x,y,z);
max_mask = (180/pi)*max(elec(:,1))-3;
welec = 1;
icount = 1;
[k,kinv,a,ainv,e] = k_and_e(welec,x,y,z);
for i = 1:image_rows
	for j = 1:image_columns
	       	w = data(icount,:);
	       	v = zeros(1,size(good_chan,2));
		v = w(good_chan);
		icount = icount+1;
		[p,q,error_check]= mateqs(welec,x,y,z,v,k,kinv,a,ainv,e);
		image_3d= interp_3d(welec,x,y,z,xs,ys,zs,p,q);	
		imagemat= mask(xs,ys,zs,image_3d,max_mask); 
		image_binary(plot_res*(i-1)+1:plot_res*(i-1)+plot_res,plot_res*(j-1)+1:plot_res*j)= imagemat;
end;
end;
if (nargin == 6)
	if ~isempty(out_root)
		out_fname = [out_root '.img'];
		fid = fopen(out_fname,'wb');
		if fid < 0
			error('couldnt open output file')
		end
	else
		fid = -1;
	end
elseif nargin < 6
	fid = put_fid('wb');
	if fid < 0
		error('couldnt open output file')
	end
end;
if fid ~= -1
	size(image_binary);
	num_written = fwrite(fid, size(image_binary), 'float');
	if num_written ~= 2
			error('Failed to write header info');
	end 
	fwrite(fid,plot_res,'float');
	num_written = fwrite(fid,image_binary','float');
	if num_written ~= size(image_binary,1)*size(image_binary,2)
		error('write-out filed')
	end;
	fclose(fid);
end;
