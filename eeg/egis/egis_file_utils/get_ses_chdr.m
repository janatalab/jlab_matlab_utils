function chdr = get_ses_chdr(fname)
% ses_chdr = get_ses_chdr(fname);
%
% Returns session cell headers from EGIS file specified in fname
%
% This function provides a wrapper around the file I/O that necessarily
% accompanies lower-level routines such as rd_egis_hdr_v.
%

% 10/22/01 PJ

byte_sex = 'ieee-be';

fid = fopen(fname, 'rb', byte_sex);

if fid == -1, error(sprintf('Error opening file %s', fname)), end

[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(fid);

fclose(fid);

return