function status=laplacian_ses(bad_chan,infname)
%function status=laplacian_ses(bad_chan,infname)
%calculate surface laplacian of a session file.  
%
%parameters
%	infname = input filename
%	outfname = output filename.  
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
	srcfid = fopen(infname, 'r');
	if srcfid == -1
		error(['Could not open input file ' infname '.']);
	end
	outfname = [infname '_lap'];
	destfid = fopen(outfname, 'wb');
	if destfid == -1
		temp =fclose(srcfid);
		error(['Could not open output file ' outfname '.']);
	end
%Otherwise run interactive mode
else
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*', 'Open Session File:');
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
ses_hdr_offsets_v;
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
[xelec,yelec,zelec] = electrodes(fhdr(NChan)+1);
ses_mask = artifact_edit(infname, fhdr, chdr,0);

%Begin looping through cells
for c=1:fhdr(NCells)
	disp(['Processing cell: ' int2str(c)]);

	%Loop through observations
	for t=1:chdr(c,NObs)
		ses_mask_cell_offset = sum(chdr(1:c, NTrials)) - chdr(c,NTrials);
		ref_trialdata = zeros(chdr(c,NPoints), fhdr(NChan));		
		if sum(ses_mask(ses_mask_cell_offset + t,:)) > min_good_chan
			trialdata = rd_onetr_allch(srcfid, coff(c), t, fhdr(NChan), chdr(c, NPoints));
			avref_trialdata= average_reference(trialdata, ses_mask(ses_mask_cell_offset + t,:));
			lap_trialdata = zeros(size(avref_trialdata,1),size(avref_trialdata,2));
			
			good_chan = find(ses_mask(ses_mask_cell_offset + t,:));
			lap_trialdata = laplacian_trial(trialdata,good_chan,xelec,yelec,zelec)
			ref_trialdata = rereference(lap_trialdata);
		end;
		fwrite(destfid, ref_trialdata', 'int16');
	end 	%for t=1:chdr(c,NObs)
end 	%for c=1:fhdr(NCells)

disp('Finished running laplacian_ses. Have a nice day');
temp=fclose(srcfid);
temp=fclose(destfid);
status = 1;
