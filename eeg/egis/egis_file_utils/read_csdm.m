function [csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdm_fname, cell, obs)

% function [csdm_data, csdm_fhdr, csdm_chdr] = read_csdm(csdm_fname, cell, obs)

size_data_type = 4;

if nargin < 3
  error('Need 3 input arguments');
end

if length(cell) > 1
  error('Only one cell can be specified at a time');
end

if length(obs) > 1
  error('Only one observation can be specified at a time');
end

% Open the csdm file

csdm_fid = fopen(csdm_fname, 'rb');

% Read the csdm file header

ave_hdr_offsets_v;

[csdm_fhdr,csdm_chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_csdm_hdr_v(csdm_fid);

[cell_data_offsets]=get_csdm_cell_off(csdm_fhdr, csdm_chdr);

% Seek to desired place in file

foffset = cell_data_offsets(cell) + (obs-1)*csdm_fhdr(NChan)*csdm_chdr(cell,NPoints)*size_data_type;

fseek(csdm_fid, foffset, 'bof');

real_csdm_data = fread(csdm_fid, [csdm_fhdr(NChan) csdm_chdr(cell, NPoints)/2], 'float');

imag_csdm_data = fread(csdm_fid, [csdm_fhdr(NChan) csdm_chdr(cell, NPoints)/2], 'float');

csdm_data = real_csdm_data' + imag_csdm_data'*i;

fclose(csdm_fid);
