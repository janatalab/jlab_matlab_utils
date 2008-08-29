function status=make_cortical(sigma_v,sigma_m,bad_chan,bad_obs,infname,outfname)
%function status=make_cortical(sigma_v,sigma_m,bad_chan,bad_obs,infname,outfname)
%
%For now just estimates cortical surface dipoles which is the same thing as a 
%Laplacian.  Will eventually integrate evidence etc. once i know how
%to use it
%
%This function accepts either 2, 4 or 6 input arguments. If the
%last two arguments are omitted the function is run in 
%interactive mode, and file names are obtained through the std
%dialog boxes. 
%
%sigma_v and sigma_m are variance (noise) estimates on the data and the 
%dipole moments respectively
%
%Uses: ave_hdr_offsets_v, rd_egis_hdr_v, zeromeans, rd_onetr_allch

%Modification History
%	9/95 Created by Jerry and Ramesh 
%	3/22/96 bad_obs filtering by BCR

status = -1;

%Check number of input arguments
if ~((nargin == 6)|(nargin == 2)|(nargin == 4))
	error('Number of input arguments must be 2, 4 or 6.');
end
if nargin == 2
	bad_chan = [];bad_obs=[];
end;
if sigma_m == 0
	sigma_m = 1;
end;
%Initialize fids
srcfid = -1;
destfid = -1;
%First try batch mode
if nargin == 6
	srcfid = fopen(infname, 'r');
	if srcfid == -1
		error(['Could not open input file' infname '.']);
	end
	destfid = fopen(outfname, 'wb');
	if destfid == -1
		temp =fclose(srcfid);
		error(['Could not open output file ' outfname '.']);
	end
%Otherwise run interactive mode
else
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*.ave*', 'Open Average File:');
	end
	while destfid == -1
		[destfid, outfname]=put_fid('wb','new.ave_crt','Save New CorticalAverage File As:');
	end
end

%Call EGIS hdr index script
ave_hdr_offsets_v;
%Read in the EGIS header
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);
%Calculate nyquist frequency and compare to lowpass parameter
% Copy source header to destination file
frewind(srcfid);
[temp, countin] = fread(srcfid, fhdr(LHeader), 'char');
countout = fwrite(destfid, temp, 'char');
if countin ~= countout
  error('Error copying header information');
end
if (nargin > 2)&(bad_chan ~= [])
	if size(bad_chan,1) ~= chdr(1,NObs)
		error('mismatch between bad_chan and NObs')
	end;
end;
[xelec,yelec,zelec] = electrodes(fhdr(NChan));
sources = xyz2tp(xelec,yelec,zelec);
A = transfer_matrix(500,0.2,80,8.0,8.2,8.7,9.2,7.8,sources,sources);

%Begin looping through cells
for c=1:fhdr(NCells)
	disp(['Processing cell: ' int2str(c)]);
	%Loop through subject averages
	for t=1:chdr(c,NObs)
		if ~any(t==bad_obs)
			trialdata = rd_onetr_allch(srcfid, coff(c), t, fhdr(NChan), chdr(c, NPoints));
			trial_mask = ones(1,fhdr(NChan));
			if (nargin > 2)&(bad_chan ~= [])
				trial_mask(1,bad_chan(t,find(bad_chan(t,:)))) = zeros(1,size(find(bad_chan(t,:)),2));
			end;
			good_chan = find(trial_mask);
			good_trial = zeros(size(trialdata,1),size(good_chan,2));
			good_trial = trialdata(:,good_chan);
			good_A = zeros(size(good_chan,2),size(good_chan,2));
			good_A = A(good_chan,good_chan);
			
			
			cortical_trial = zeros(size(trialdata,1),size(trialdata,2));
			[m,deviations] = bayes_dipole_trial(good_A,good_trial',sigma_v,sigma_m);
			cortical_trial(:,good_chan) = m';
			cortical_trial = cortical_trial*10;
			fwrite(destfid,cortical_trial', 'int16');
		else
			cortical_trial=zeros(chdr(c, NPoints),fhdr(NChan));
			fwrite(destfid,cortical_trial', 'int16');
		end
	end 	%for t=1:chdr(c,NObs)
end 	%for c=1:fhdr(NCells)

disp('Welcome to the cortical surface.  Transfer function fee of $5 applies');
temp=fclose(srcfid);
temp=fclose(destfid);
status = 1;
