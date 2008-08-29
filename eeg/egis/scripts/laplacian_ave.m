function status=laplacian_ave(bad_chan,infname)
%function status=laplacian_ave(bad_chan,infname) 
%calculate surface laplacian of an average ERP file.  Do not use this 
%on power or coherence files.  
%
%parameters
%	bad_chan = bad_channels (optional).  If bad channels are specified
%	there must be one row for every observation in the egis file
%	infname = input filename (optional)

%Modification History
%	2/96 Created by Ramesh Srinivasan, master of the universe
%
%comments on modification history
%
status = -1;

%Check number of input arguments
if ~((nargin == 1)|(nargin == 2)|(nargin == 0))
	error('Number of input arguments must be either 1 2 or 0 moron.');
end
%Initialize fids
srcfid = -1;
destfid = -1;
if nargin == 0
	bad_chan = [];
end;
%First try batch mode
if nargin == 2
	srcfid = fopen(infname, 'r','b');
	if srcfid == -1
		error(['Could not open input file ' infname '.']);
	end
	outfname = [infname '_lap'];
	destfid = fopen(outfname, 'wb','b');
	if destfid == -1
		temp =fclose(srcfid);
		error(['Could not open output file ' outfname '.']);
	end
%Otherwise run interactive mode
else
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*', 'Open Average File:');
	end
	outfname = [infname '_lap'];
	destfid = fopen(outfname, 'wb');
	if destfid == -1
		temp =fclose(srcfid);
		error(['Could not open output file ' outfname '.']);
	end
	end
end

%Call EGIS hdr index script
ave_hdr_offsets_v;
%Read in the EGIS header
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);


% Copy source header to destination file
frewind(srcfid);
[temp, countin] = fread(srcfid, fhdr(LHeader), 'char');
countout = fwrite(destfid, temp, 'char');
if countin ~= countout
  error('Error copying header information');
end
% CREATE CHANNEL MASK
good_chan_mask = ones(chdr(1,NObs),fhdr(NChan));


if bad_chan ~= []
  % Check to see if number of rows in obs and bad_chans jives
  nrows_bad = size(bad_chan,1);
  if chdr(1,NObs) ~= nrows_bad
    disp(['WARNING: Number of observations (' int2str(nobs) ') doesn''t match number of rows in bad_chan (' int2str(nrows_bad) ')']);
    if nrows_bad == 1
      disp('Assuming bad_chan applies to all observations');
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

%Begin looping through cells
for c=1:fhdr(NCells)
	disp(['Processing cell: ' int2str(c)]);

	%Loop through subject averages
	for t=1:chdr(c,NObs)
		trialdata = rd_onetr_allch(srcfid, coff(c), t, fhdr(NChan), chdr(c, NPoints));
		lap_trialdata = zeros(size(trialdata,1),size(trialdata,2));
		good_chan = find(good_chan_mask(t,:));
		lap_trialdata = laplacian_trial(trialdata,good_chan,xelec,yelec,zelec);
		fwrite(destfid, lap_trialdata', 'int16');
		
		
	end 	%for t=1:chdr(c,NObs)
end 	%for c=1:fhdr(NCells)

disp('Finished running laplacian_ave. Have a nice day');
temp=fclose(srcfid);
temp=fclose(destfid);
status = 1;