function [ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext] = csdmhdr2avehdr(csdmfiles, obs)
%  function [ave_fhdr,ave_chdr,ave_ename,ave_cnames,ave_fcom,ave_ftext] = csdmhdr2avehdr(csdmfiles, obs)
%
%  Moves relevant header information from csdm header to average header
%
%  csdmfiles -- Matrix of csdm filenames
%  obs -- Observations to extract from csdm files
%
%  Note:  	In the event that there is more than one csdmfile and length(obs) > 1,
%			multiple files will be written out, with one file/obs.

%  Modification history:
%  01/26/96 PJ	Started work on function

%
% read in a csdm header in order to make an output file header
%
ave_hdr_offsets_v;

if nargin < 2
  obs = [];
end

csdm_fid = fopen(csdmfiles(1,:),'r');
if csdm_fid == -1
  error(['Trouble opening file -- ' csdmfiles(1,:)])
end

[csdm_fhdr,csdm_chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_csdm_hdr_v(csdm_fid);
fclose(csdm_fid);

if obs == []
  obs = csdm_chdr(1,NObs);
end

if any(obs > csdm_chdr(1,NObs))
	error('observation number out of range')
end;

ave_fhdr = csdm_fhdr;
if csdm_fhdr(NChan) == 8385
	ave_fhdr(NChan) = 129;
elseif csdm_fhdr(NChan) == 2145
	ave_fhdr(NChan) = 65;
elseif csdm_fhdr(NChan) == 33153;
	ave_fhdr(NChan) = 257;
else
	error('Nchan is unacceptable in this csdm file')
end;

% Copy average file subject specifics array from each input file
ave_chdr = concat_csdm_chdr(csdmfiles, obs);

ave_ename = ename; ave_fcom = fcom; ave_ftext = ftext;
ave_cnames = cnames;

ave_fhdr(LastDone) = size(csdmfiles,1);

% Write the correct number of observations into new header
% depending on  

for c = 1:csdm_fhdr(NCells)
  if length(obs) > 1 & size(csdmfiles,1) == 1
    ave_chdr(c,NObs) = length(obs);
  else
    ave_chdr(c,NObs) = size(csdmfiles,1);
  end
end
