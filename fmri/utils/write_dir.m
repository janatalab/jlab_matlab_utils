function write_dir(V,dirname,fstub)
% write_dir(V,dirname,fstub);
%
% Writes image volume time series data from V into directory dirname
% fstub -- root of individual filenames.  Final file names are
%
%    fstub_i0001.img, fstub_i0002.img, etc.
%

% 5/18/00 PJ

% Check to see if dirname exists

if ~exist(dirname,'dir')
  disp(sprintf('Creating directory: %s', dirname))
  unix(['mkdir ' dirname]);
end

dims = size(V);
nvol = dims(4);

disp(sprintf('Writing %d volumes ...', nvol))

for ivol = 1:nvol
  outfname = sprintf('%s/%s_i%04d.img', dirname, fstub, ivol);
  
  fid = fopen(outfname, 'wb', 'l');
  if fid == -1
    error(sprintf('Failed to open file %s for writing', outfname))
  end
  
  count = fwrite(fid,V(:,:,:,ivol),'integer*2');
  if count ~= prod(dims(1:3))
    error(sprintf('Failed to write all data for file %s', outfname))
  end
    
  fclose(fid);
end % for ivol=
