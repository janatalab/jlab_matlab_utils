function [fhdr, chdr] = get_csdm_hdr(fname)
%function [fhdr, chdr] = get_csdm_hdr(fname);
%
%

csdm_fid = fopen(fname,'r');
if csdm_fid == -1
  error(['Trouble opening file -- ' fname])
end
[fhdr,chdr]=rd_csdm_hdr_v(csdm_fid);
fclose(csdm_fid);
