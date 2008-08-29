function new_chdr = concat_csdm_chdr(csdmfiles, obs)
% function new_chdr = concat_csdm_chdr(csdmfiles, obs)
%
% Concatenates several average file subject specifics array
% into a single chdr
% NOTE: THIS ROUTINE ASSUMES HOMOGENEOUS CELL HEADERS

% Modification history:
% 01/26/96  PJ Created routine
%

csdm_hdr_offsets_v;

csdm_fid = fopen(csdmfiles(1,:),'r');
[csdm_fhdr,csdm_chdr]=rd_csdm_hdr_v(csdm_fid);
fclose(csdm_fid);

num_files = size(csdmfiles,1);

if nargin == 2 & num_files > 1
  if length(obs) > 1
    disp('Warning: Multiple files, multiple observations')
    disp('To have proper headers, either specify one observation, or use all, dammit!')
  end
else
  obs = 1:csdm_chdr(1,NObs);
end

num_obs = length(obs);

% Copy top part of chdr from 1st file
new_chdr = csdm_chdr(:,1:5);

for f = 1:size(csdmfiles,1)
  csdm_fid = fopen(csdmfiles(f,:),'r');
  [csdm_fhdr,csdm_chdr]=rd_csdm_hdr_v(csdm_fid);
  fclose(csdm_fid);

  last_col = size(new_chdr,2);
  num_new_col = num_obs * csdm_chdr(1,LSpec)/2;
  if num_obs > 1
    new_chdr(:,last_col+1:last_col+num_new_col) = csdm_chdr(:,6:5+num_new_col);
  else
    spec_start = 6+(obs-1)*num_new_col;
    spec_stop = 5 + (obs)*num_new_col;
    new_chdr(:,last_col+1:last_col+num_new_col) = csdm_chdr(:,spec_start:spec_stop);
  end
end

return
