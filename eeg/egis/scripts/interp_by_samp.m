function status = interp_by_samp(cells,obs,samps,plot_res,interp_type,bad_chan,output_root,header_flag,infname,stack)
%function status = interp_by_samp(cells,obs,samps,plot_res,interp_type,bad_chan,output_root,header_flag,infname, stack)
%
%cells = rows by columns array of cell mapping to page
%obs = rows by columns array of obs mapping to page
%samps = vector of samples to interpolate
%plot_res = grid size
%interp_type - 'sp' = spherical harmonic = '3d' for 3-d splines
% bad_chan = array or matrix of bad_channels 
% output_root = root name for output image
% header_Flag = flag top determine swhether the images get a header, 'h' for header 's' for no header defaults of 'h'
% stack = if this is true, interpolations will be placed into 1 file

%Modification History
%	2/95 Created by Ramesh Srinivasan
%		  included features from PJ's update to make_interp_by_samp which is now defunct
%  02/23/96  PJ  nobs variable is set equal to number of rows in obs.  There is
%                still a fundamental problem with the way obs and bad_chan
%                are handled.
%
%  01/21/97 PJ Added option to put all interpolations into a single file
%
%  04/26/97 PJ Updated to handle more efficient routines in /spline
%              Handling of more than one obs is almost certainly defunct.
%
%  11/23/97 PJ Fixed masking problem under Matlab 5.0

%Check number of input arguments
if ~((nargin ==9)|(nargin == 5)|(nargin == 6)|(nargin == 8)|(nargin == 10))
	error('Number of input arguments must be either 5 6 8 or 9 or 10.');
end
%Initialize fids
srcfid = -1;

if ~((size(obs,1) == size(cells,1))&(size(obs,2) == size(cells,2)))
	error('page layout is faulty')
end;
if nargin == 5
	bad_chan = [];
	header_flag = 'h';
	output_root = 'interpolation';
end;
if isempty(output_root)
	output_root = [infname '_interp_'];
end;
if ~((header_flag == 'h')|(header_flag == 's'))
	header_flag = 'h';
end;
if isempty(samps),  samps = 1:chdr(1:NPoints), end;

if nargin < 10
  stack = 0;
end

%First try batch mode
if nargin == 9 | nargin == 10
	srcfid = fopen(infname, 'r','ieee-be');
	if srcfid == -1
		error(['Could not open input file' infname '.']);
	end

%Otherwise run interactive mode
else
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*.*', 'Open File:');
	end
end

%Call EGIS hdr index script
ave_hdr_offsets_v;

%Read in the EGIS header
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);

% CREATE CHANNEL MASK
good_chan_mask = ones(chdr(1,NObs),fhdr(NChan));

nobs = size(obs,1);

if ~isempty(bad_chan)
  % Check to see if number of rows in obs and bad_chans jives
  nrows_bad = size(bad_chan,1);
  if nobs ~= nrows_bad
    disp(['WARNING: Number of observations (' int2str(nobs) ') doesn''t match number of rows in bad_chan (' int2str(nrows_bad) ')']);
    if nrows_bad == 1
      disp('Replicating bad_chan vector to match nobs');
      bad_chan = ones(size(good_chan_mask,1))*bad_chan;
    else
      error('Drastic mismatch between nobs and bad_chan')
    end
  end

  for isub = 1:size(good_chan_mask,1)
    bchan = find(bad_chan(isub,:));
    good_chan_mask(isub,bad_chan(isub,bchan)) = zeros(1,size(bchan,2));
  end;
end;

[xelec,yelec,zelec] = electrodes(fhdr(NChan));
if interp_type == 'sp'
	if fhdr(NChan) == 129
		nmax = 8;
	else
		error('something is rotten')
	end;
end;

good_chan = ones(1,fhdr(NChan));
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

% WORK THRHOUGH EACH SAMPLE

for isamp = 1:size(samps,2)
  disp(['Working on sample: ' int2str(samps(isamp))]);
  image_binary = zeros(plot_res*size(cells,1),plot_res*size(cells,2));
  for icell=1:size(cells,1)
	    for t=1:size(cells,2)
	    c = cells(icell,t);
%	    disp(['  Cell: ' int2str(c)]);
%	    disp(['    Obs: ' int2str(obs(icell,t))]);
   	    trialdata = rd_onetr_allch(srcfid, coff(c), obs(icell,t), fhdr(NChan), chdr(c, NPoints));
        w = trialdata(samps(isamp),:);
        v = w(good_chan);

		if interp_type == 'sp'
	        meanv = mean(v);
    	    v = v - meanv*ones(1,size(v,2));
    		  [sp_coff,error_check, sp_mat] = bayes_sphcoff(nmax,x,y,z,v);
    		  imagemat = image_2D(sp_coff,nmax,plot_res,120,meanv);
		elseif interp_type == '3d'
	  		[p,q,error_check]= mateqs(welec,x,y,z,v,k,kinv,a,ainv,e);
			ruu = sqrt(x(1).^2+y(1).^2+z(1).^2);
			[xs,ys,zs] = polgrid(plot_res,ruu);	
			image_3d= interp_3d(welec,x,y,z,xs,ys,zs,p,q);	
			imagemat= mask(xs,ys,zs,image_3d,120); 
		else
			error('havent dealt with spherical splines yet')
		end;			
    	image_binary(plot_res*(icell-1)+1:plot_res*(icell-1)+plot_res,plot_res*(t-1)+1:plot_res*(t-1)+plot_res)=imagemat;
 		
    end;  %for t=1:size(obs,2)
  end 	%for c=1:size(cells,2)

  if ~stack  
    outfname = [output_root '_s' int2str(samps(isamp))];
    outfid = fopen(outfname,'wb','ieee-be');
  elseif isamp == 1
    outfname = [output_root '_s' int2str(samps(isamp)) '_s' int2str(samps(length(samps)))]
    outfid = fopen(outfname,'wb','ieee-be');
  end
  
  if isamp == 1 | ~stack
    if header_flag == 'h'    
	  num_written = fwrite(outfid, size(image_binary), 'float');
  	  if num_written ~= 2
    	error('Failed to write header info');
  	  end

  	fwrite(outfid, plot_res, 'float');
    end;
  end
  	num_written = fwrite(outfid,image_binary','float');
  	if num_written ~= size(image_binary,1)*size(image_binary,2)
    	error('write-out filed')
  	end;
  if ~stack
    fclose(outfid);
  end
end 	%for c=1:size(samps,2)

disp('Finished running make_interp.');
temp=fclose('all');
status = 1;
