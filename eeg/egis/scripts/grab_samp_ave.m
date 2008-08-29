function trialdata=grab_samp_ave(cell,obs,infname)
%function status=make_ave_bf(prestimsamps,prestimsamps2,lowpass,highpass,infname,outfname)
%
%Filter and baseline correct data in an average EGIS file.
%
%This function accepts either 4 or 6 input arguments. If the
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
status = -1;

%Check number of input arguments
if ~((nargin == 3))
	error('Number of input arguments must be 2');
end
%Initialize fids
srcfid = -1;
srcfid = fopen(infname,'r')
%Call EGIS hdr index script
ave_hdr_offsets_v;
%Read in the EGIS header
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid);
%Calculate nyquist frequency and compare to lowpass parameter
%Begin looping through cells
c= cell;
t = obs;
trialdata = rd_onetr_allch(srcfid, coff(c), t, fhdr(NChan), chdr(c, NPoints));
temp=fclose(srcfid);
temp=fclose(destfid);
status = 1;
