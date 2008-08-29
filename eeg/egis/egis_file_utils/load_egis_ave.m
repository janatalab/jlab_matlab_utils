function ave = load_egis_ave(fname)
% ave = load_egis_ave(fname);
%
% Loads entire EGIS .ave or .gave file
%


% 01/24/00 PJ
byte_sex = 'ieee-be';

ave_hdr_offsets_v;
scaling_factor = 500;

fid = fopen(fname, 'rb', byte_sex);
[d.fhdr,d.chdr,d.ename,d.czeros,d.cgains,d.cnames,d.fcom,d.ftext, coff]=rd_egis_hdr_v(fid);

d.data = zeros(d.chdr(1,NPoints),d.fhdr(NChan),size(d.chdr,1));

for icell = 1:size(d.chdr,1)
  d.data(:,:,icell) = rd_onetr_allch(fid, coff(icell), 1, d.fhdr(NChan), d.chdr(icell,NPoints))/scaling_factor;
end

fclose(fid);

ave = d;
