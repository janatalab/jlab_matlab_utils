function [celldata] = rd_onecell_allchan(fid, fhdr, chdr, cellnum)
%function [celldata] = rd_onecell_allchan(fid, fhdr, chdr, cellnum)
%
%Read in all the data for a given cell
%

if fhdr(2) == -1
  ave_hdr_offsets_v;
else
  ses_hdr_offsets_v;
end

celldata = zeros(chdr(cellnum, NObs)*chdr(cellnum, NPoints), fhdr(NChan));

cell_data_offsets = get_cell_offsets(fhdr, chdr);

for t = 1:chdr(cellnum, NObs)
	celldata((t-1)*chdr(cellnum, NPoints)+1:t*chdr(cellnum, NPoints),:) = rd_onetr_allch(fid, cell_data_offsets(cellnum), t, fhdr(NChan), chdr(cellnum,NPoints));
end
