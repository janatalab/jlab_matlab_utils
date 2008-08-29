function img_stack = interp_no_egis(data,bad_chan,interp_type,plot_res,outfile,header_flag,stack)
% function img_stack = interp_no_egis(data,bad_chan,interp_type,plot_res,outfile,header_flag, stack);
%
% Calls interpolation routines but doesn't rely on EGIS file for
% the data
%
% Adapted from existing interpolation routines written by RS
%
% img_stack is size plot_res,plot_res,size(data,1)

% Modification history:
% 11/30/96 PJ Started script
% 06/05/00 PJ Modified further so that if output file arguments aren't passed,
%             routine simply returns data in img_stack.
%

write_images = 1;

if nargin < 1
  disp(['Not enough input arguments to interp_no_egis'])
  return
elseif nargin < 2
  bad_chan = [];
  interp_type = '3d';
  plot_res = 40;
  outfile = 'interp.out';
  header_flag = 'h';
  stack = 1;
elseif nargin < 3
  interp_type = '3d';
  plot_res = 40;
  outfile = 'interp.out';
  header_flag = 'h';
  stack = 1;
elseif nargin < 4
  plot_res = 40;
  outfile = 'interp.out';
  header_flag = 'h';
  stack = 1;
elseif nargin < 5
  write_images = 0;
elseif nargin < 6
  header_flag = 'h';
  stack = 1;
elseif nargin < 7
  stack = 1;
end

if write_images & ~stack
  disp(['Only supporting stacked output at this time'])
  return
end

if any(size(data) == 1)
  nsamps = 1;
else
  nsamps = size(data,1);
end

good_chan = ones(1,size(data,2));
if ~isempty(bad_chan)
  good_chan(bad_chan) = zeros(size(bad_chan));
end

good_chan = find(good_chan);

[xelec,yelec,zelec] = electrodes(129);

  v = zeros(1,size(good_chan,2));
  x = zeros(1,size(good_chan,2));
  y = zeros(1,size(good_chan,2));
  z = zeros(1,size(good_chan,2));

  x = xelec(good_chan);
  y = yelec(good_chan);
  z = zelec(good_chan);

welec = 1;

[k, kinv, a, ainv, e] = k_and_e(welec,x,y,z);

img_stack = zeros(plot_res,plot_res,nsamps);

for isamp = 1:nsamps

  disp(['  samp: ' int2str(isamp)])

  w = data(isamp,:);
  v = w(good_chan);

  [p,q]= mateqs(welec,x,y,z,v,k,kinv,a,ainv,e);
%  disp('plot is being made as a 2d image')
  ruu = sqrt(x(1).^2+y(1).^2+z(1).^2);
  [xs,ys,zs] = polgrid(plot_res,ruu);	
  image_3d= interp_3d(welec,x,y,z,xs,ys,zs,p,q);	
  imagemat= mask(xs,ys,zs,image_3d,120); 

  img_stack(:,:,isamp) = imagemat;
  
  if write_images
    if isamp == 1
      outfname = [outfile '_s' int2str(isamp) '_s' int2str(nsamps)];
      outfid = fopen(outfname,'wb','ieee-be');
      
      if header_flag == 'h'    
	num_written = fwrite(outfid, size(imagemat), 'float');
	if num_written ~= 2, error('Failed to write header info'), end
	fwrite(outfid, plot_res, 'float');
      end
    end
    
    num_written = fwrite(outfid,imagemat','float');
    if num_written ~= size(imagemat,1)*size(imagemat,2), error('write-out failed'), end
  end
  
end % for isamp=

if write_images
  fclose(outfid);
end

return





