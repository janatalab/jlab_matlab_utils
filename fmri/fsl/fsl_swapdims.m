function status = fsl_swapdims(infname,swaplist,outfname,verbose)
% status = fsl_swapdims(infname,swaplist,outfname)
%
% Swaps dimensions in file infname according to the list given in swaplist. 
% The function uses fslswapdim which is part of the FSL distribution.
%
% Example: 
% fsl_swapdims(orig_hires,{'z -x y','-x y z'}, new_hires)
%
% This would call fslswapdim twice, the first time placing 1) the z-dimension of
% the original image in to the x-dimension of a temporary image, 2) a flipped
% version of the x-dimension of the original image into the y-dimension of a
% temporary image, and 3) the y-dimension of the original image into the
% z-dimension of the temporary image.  The second swap specified by, '-x y z',
% would create a final image given by outfname.  If outfname is not specified,
% the original image is replaced!
% 
% 02/16/06 Petr Janata
% 03/23/08 Fred Barrett - added additional handling of nii.gz files
% generated by avwswapdim, when executing avwchfiletype
% 08/17/08 FB - seeded the random number generator within this function
% 2009.03.11 FB - changed from avw* to fsl* utilities, with move to fsl v4.1

status = 0;

if nargin < 3
  fprintf('Will replace original image!!\n');
  outfname = infname;
end

if nargin < 4
  verbose = 0;;
end

nswap = length(swaplist);

% seed random number generator
rand('state',sum(100*clock));

% Execute the series of swaps
tmp{1} = infname;
for iswap = 1:nswap
  tmp{iswap+1} = sprintf('/tmp/tmp_%06d', fix(rand*999999));
  fsl_str = sprintf('fslswapdim %s %s %s', tmp{iswap}, swaplist{iswap}, tmp{iswap+1});
  if verbose, fprintf('%s\n', fsl_str); end
  status = unix(fsl_str);   
  if status, return, end

  % Generate the full name of the temporary output file
  unix_str = sprintf('ls %s.*', tmp{iswap+1});
  if verbose, fprintf('%s\n', unix_str); end
  [status, tmpfname] = unix(unix_str);
  tmp{iswap+1} = deblank(tmpfname);  % remove any newline characters or spaces

end

% Move the last file to the desired destination and clean up any intermediate
% files

% Check to see if the file type of intermediate file and destination match
[destpath,deststub,destfmt] = fileparts(outfname);
[tmppath,tmpstub,tmpfmt] = fileparts(tmp{end});

% Convert if necessary
if isempty(strmatch(destfmt,tmpfmt,'exact'))
  switch destfmt
    case '.nii'
      dest_type = 'NIFTI';
    case {'.img','.hdr'}
      dest_type = 'NIFTI_PAIR';
  end

  fsl_str = sprintf('fslchfiletype %s %s', dest_type, tmp{end});
  if verbose, fprintf('%s\n', fsl_str); end
  status = unix(fsl_str);
  if status, return, end
  
  if (strmatch('.gz',tmpfmt,'exact'))
      [tmppath0,tmpstub,tmpfmt0] = fileparts(tmpstub);
  end
  
  tmp{end} = fullfile(tmppath, [tmpstub destfmt]);
end

% Move the temporary output file to the destination output file
unix_str = sprintf('mv -f %s %s', tmp{end},outfname);
if verbose, fprintf('%s\n', unix_str); end
status = unix(unix_str);
if status, error, end
 
% If we are dealing with an image pair, make sure we move the .hdr also
if any(strmatch(destfmt,{'.img'},'exact'))
  unix_str = sprintf('mv -f %s %s', fullfile(tmppath, [tmpstub '.hdr']), fullfile(destpath,[deststub '.hdr']));
  if verbose, fprintf('%s\n', unix_str); end
  status = unix(unix_str);
  if status, error, end
end

% Delete the intermediate files, but not the original file!
for iswap = 1:nswap-1
  unix_str = sprintf('rm %s', tmp{iswap+1});
  if verbose, fprintf('%s\n', unix_str); end
  status = unix(unix_str);
  if status, return, end
end

% Now we have to delete what might be erroneous orientation information
fsl_str = sprintf('fslorient -deleteorient %s', outfname);
if verbose, fprintf('%s\n', fsl_str); end
status = unix(fsl_str);   
if status, return, end
