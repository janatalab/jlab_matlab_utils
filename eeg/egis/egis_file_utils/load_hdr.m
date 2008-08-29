function [fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext,coff]=load_hdr(fname)
%[fhdr,chdr,...,coff]=load_hdr(fname)
%
%This function opens an EGIS file, reads its header,
%returns the values from rd_egis_hdr_v, then closes
%then file.
%
%See also: rd_egis_hdr_v.

fid=fopen(fname);
if fid == -1, error('Fuck up opening input file.'), end
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext,coff]=rd_egis_hdr_v(fid);
status=fclose(fid);

