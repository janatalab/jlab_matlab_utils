function [cell_data_offsets]=get_csdm_cell_off(fhdr, chdr)
%function [cell_data_offsets]=get_csdm_cell_off(fhdr, chdr)
%
%Returns the index to the beginning of each cell's data
%
% Modification history:
% 1/10/96 PJ Adapted from get_cell_offsets, but now multiplying by
%            size_data_value to reflect fact that we are writing floats

ave_hdr_offsets_v;

size_data_value = 4;  % 4=float

cell_data_offsets = zeros(1,fhdr(NCells));
cum_cell_sizes = 0;

for c = 1:fhdr(NCells)
	cell_data_offsets(c) = fhdr(LHeader) + cum_cell_sizes;
	cum_cell_sizes = cum_cell_sizes + chdr(c,NObs)*chdr(c,NPoints)*fhdr(NChan)*size_data_value;
end

