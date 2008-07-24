function [img,dims,scales,bpp,endian] = read_avw(fname)
% [img, dims,scales,bpp,endian] = READ_AVW(fname)
%
%  Read in an Analyze or nifti file into either a 3D or 4D
%  array (depending on the header information)
%  fname is the filename (must be inside single quotes)
%  Note: automatically detects - unsigned char, short, long, float
%         double and complex formats
%  Extracts the 4 dimensions (dims), 
%  4 scales (scales) and bytes per pixel (bpp) for voxels 
%  contained in the Analyze or nifti header file (fname)
%  Also returns endian = 'l' for little-endian or 'b' for big-endian
%
%  See also: SAVE_AVW

% remove extension if it exists
if ( (length(findstr(fname,'.hdr'))>0) | ...
        (length(findstr(fname,'.img')>0)) | ...
        (length(findstr(fname,'.nii')>0)) ),
  fname=fname(1:(length(fname)-4));
end
if ( (length(findstr(fname,'.hdr.gz')>0)) | ...
        (length(findstr(fname,'.img.gz')>0)) | ...
        (length(findstr(fname,'.nii.gz')>0)) ),
  fname=fname(1:(length(fname)-7));
end

%% convert to uncompressed nifti pair (using FSL)
tmpname = tempname;
tmp=sprintf('! $FSLDIR/etc/matlab/convertTOniftipair %s %s \n',fname,tmpname);
eval(tmp);

  [dims,scales,bpp,endian,datatype]= read_avw_hdr(tmpname);
  if (datatype==32),
    % complex type
    img=read_avw_complex(tmpname);
  else
    img=read_avw_img(tmpname);
  end
  
tmp=sprintf('! /bin/rm %s.hdr %s.img \n',tmpname,tmpname);
eval(tmp);

