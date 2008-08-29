function status=make_ave_bf(prestimsamps1,prestimsamps2,lowpass,highpass,lowfiltord,highfiltord,infname,outfname)
%function status=make_ave_bf(prestimsamps,prestimsamps2,lowpass,highpass,lowfiltord,highfiltord,infname,outfname)
%
%Filter and baseline correct data in an average EGIS file.
%
%This function accepts either 4, 6 or 8 input arguments. If the
%last two arguments are omitted the function is run in 
%interactive mode, and file names are obtained through the std
%dialog boxes. 
%
%The parameter prestimsamps1 & 2 are used to baseline correct through
%a call to zeromeans (c.f.). The parameters lowpass and highpass
%are used to construct the butterworth filter. If highpass is zero
%than filtering occurs through a single lowpass filtering of order
%20. Otherwise, data is bandpass filtered with a filter of order six,
%and then refiltered with a lowpass filter of order 20 to improve
%high-end rolloff.
%
%Uses: ave_hdr_offsets_v, rd_egis_hdr_v, zeromeans, rd_onetr_allch

%Modification History
%	9/95 Created by Ramesh Srinivasan
%	10/95 Rewrote user interface and added high pas filtering - B. Rakitin
%	10/31/95 Bug fixed in file opening branch - B. Rakitin
%   11/17/95 reordered filters in bad-pass case - RS
%	2/26/96 added pass params for filt orders - BCR

status = -1;

%Check number of input arguments
if ~((nargin == 4)|(nargin == 6)|(nargin == 8))
	error('Number of input arguments must be either 4, 6, or 8.');
end
%Initialize fids
srcfid = -1;
destfid = -1;
%First try batch mode
if nargin == 8
	srcfid = fopen(infname, 'r');
	if srcfid == -1
		error(['Could not open input file' infname '.']);
	end
	destfid = fopen(outfname, 'wb');
	if destfid == -1
		temp =fclose(srcfid);
		error(['Could not open output file ' outfname '.']);
	end
	if highfiltord==[]
		highfiltord=3;
	end
	if lowfiltord==[]
		lowfiltord=20;
	end
%Otherwise run interactive mode
else
	while srcfid == -1
		[srcfid, infname]=get_fid('r','*ave*', 'Open Average File:');
	end
	while destfid == -1
		[destfid, outfname]=put_fid('wb','new.ave_bf','Save New Average File As:');
	end
	if nargin==6
		if highfiltord==[]
			highfiltord=3;
		end
		if lowfiltord==[]
			lowfiltord=20;
		end
	else
		highfiltord=3;
		lowfiltord=20;
	end
end

%Call EGIS hdr index script
ave_hdr_offsets_v;
%Read in the EGIS header
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);
%Calculate nyquist frequency and compare to lowpass parameter
nyquist=zeros(1,fhdr(NCells));
for c=1:fhdr(NCells)
	nyquist(c) = chdr(c, SampRate)/2;
end
if any(lowpass > nyquist)
	temp=fclose(srcfid);
	temp=fclose(destfid);
	error(['The lowpass freq. exceeds nyquist freq., ' int2str(nyquist) ', in some cells.']);
end


% Copy source header to destination file
frewind(srcfid);
[temp, countin] = fread(srcfid, fhdr(LHeader), 'char');
countout = fwrite(destfid, temp, 'char');
if countin ~= countout
  error('Error copying header information');
end

%Begin looping through cells
for c=1:fhdr(NCells)
	disp(['Processing cell: ' int2str(c)]);
	%Construct the filters
	if highpass ~= 0
		%filter_order = 3;
		cutoff = [(highpass/nyquist(c)) (lowpass/nyquist(c))];
		b_bandpass=[]; a_bandpass=[];
		[b_bandpass, a_bandpass] = butter(highfiltord, cutoff);
			if c == 1
				freqz(b_bandpass,a_bandpass);
				figure
			end;
				
	end;
%	construct low-pass filter
		%filter_order = 20;
		cutoff = lowpass/nyquist(c);
		b_lowpass=[]; a_lowpass=[];
		[b_lowpass,a_lowpass] = butter(lowfiltord,cutoff);
			if c == 1
				freqz(b_lowpass,a_lowpass);
			end;

	%Loop through subject averages
	for t=1:chdr(c,NObs)
		trialdata = rd_onetr_allch(srcfid, coff(c), t, fhdr(NChan), chdr(c, NPoints));
		for ch=1:fhdr(NChan)
				filtdata(:,ch) = filtfilt(b_lowpass, a_lowpass, trialdata(:,ch));
		end;
	
		if highpass ~= 0
			for ch=1:fhdr(NChan)
				filtdata(:,ch) = filtfilt(b_bandpass, a_bandpass, filtdata(:,ch));
			end;
		end;
		bfiltdata= zeros(size(filtdata,1),size(filtdata,2));
		bfiltdata = zeromean(filtdata,prestimsamps1,prestimsamps2);
		fwrite(destfid, bfiltdata', 'int16');
	end 	%for t=1:chdr(c,NObs)
end 	%for c=1:fhdr(NCells)

disp('Finished running make_ave_bf.');
temp=fclose(srcfid);
temp=fclose(destfid);
status = 1;
