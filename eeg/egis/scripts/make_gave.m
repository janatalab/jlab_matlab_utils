function [gave_data] = make_gave(nfiles,bad_chan,bad_obs,outfname,infname, ap)
%function [gave_data] =make_gave(nfiles,bad_chan,bad_obs,outfname,infname, params)
%
%version 2.0
%Grand-average data in EGIS average files.
%
%This function accepts either 1, 2, 3, or 5 arguments. The first
%argument is the number of files to be averaged together. If you specify
%1 as nfiles, then it is assumed that there are multiple subjects per
%cell. If nfiles is greater than one, it is assumed the files are single
%subject average files. In this case the EGIS grand average file 
%header is taken from the first input filename specified through either
%of the two means availible.
%
%If either outfname or infname are to be passed to this function, then
%both must be given, and the function is fully batch processed without
%addition user input. The infname matrix must have nfile rows, and each
%row must be a string of equal length. If these two arguments are not given
%then the function runs interactively, querying for an output filename and
%nfiles input files.
%
%The argument bad_chan is a matrix with nfiles rows and no more columns 
%than there were channels of acquisition used to exclude globally bad 
%channels as determined by previous editing. If nfiles is 1 then rows 
%correspond to observations otherwise rows correspond to the file in the
%equivalent row of infname.
%
%The argument bad_obs is used only if nfile is 1. In this case bad_obs
%is an array specifying observations (i.e. subjects) to be excluded from the
%grand average for all cells.

%Modification History
%
%	written 9/17/95 p.m. RS
%	defogged 9/18/95 a.m. RS
%	bad_obs added and error checking expanded 11/27/95 p.m. BCR
%	Code fixed for multi-file case, by BCR, 6/23/96
%	Bad_Obs warning removed, 10/3/96 BCR
%
%  05/31/00 PJ Cleaned up some of the poor matlab coding.
%              Since we are dealing with average data which are relatively
%              small, load all data at once.  Doing this also enables
%              conversion to z-scores (see below).
%
%              Added option to convert average data into z-scores prior to
%              averaging.  Added params variable to input.  This is a structure
%              that is intended to contain future parameters.  In this case it
%              has a field called compute_z.
%
%  06/05/00 PJ Added detrending option

byte_order = 'ieee-be';

gave_data = [];

%Check number of input arguments
%if ~((nargin == 2)|(nargin == 3)|(nargin == 1)|(nargin == 5))
%  error('Number of input arguments must be either 1, 2, 3, or 5.');
%end

%Check for correct use of Bad_Obs
%if (nfiles > 1) & (bad_obs ~= [])
%	disp('Warning: bad_obs argument is ignored when nfiles is greater than zero.');
%end

%Check input argument dimensions
%if (nargin == 5)
%	if ((size(infname,1) ~= nfiles) | (size(outfname,1) ~= nfiles))
%		error('Indicated number of files not equal to number of rows in one or more input arguments.');
%	end
%end
%Code fixed for multi-file case, by BCR, 6/23/96

if nargin < 6
  ap = struct('compute_z',0,'detrend',0');
else
  if ~isfield(ap,'compute_z'), ap.compute_z = 0, end
  if ~isfield(ap,'detrend'), ap.detrend = 0, end
end

if (nargin >= 5)
  if (size(infname,1) ~= nfiles)
    error('Indicated number of files not equal to input parameter nfiles.');
  end
  if (size(outfname,1) ~= 1)
    error('Outfname parameter must be a single row vector.');
  end
end

%Initialize fids
destfid = -1;
%First try batch mode
if nargin >= 5
  for i = 1:size(infname,1)
    srcfid(i) = fopen(infname(i,:), 'r', byte_order);
    if srcfid(i) == -1
      error(sprintf('Could not open input file %s', infname(i,:)));
    end
  end;
  destfid = fopen(outfname, 'wb', byte_order);
  if destfid == -1
    temp =fclose(srcfid);
    error(['Could not open output file ' outfname '.']);
  end
  %Otherwise run interactive mode
elseif nargin <= 3
  if nargin < 2, bad_chan = zeros(nfiles,1);, end
  if nargin < 3, bad_obs = [];, end
  for i=1:nfiles
    [srcfid(i), infname(i,:)]=get_fid('r','*.ave*', 'Open Average File:');
  end
  while destfid == -1
    [destfid, outfname]=put_fid('wb','new.gave_bf','Save Grand Average File As:');
  end
end

nfiles = size(srcfid,2);

%
% Call EGIS hdr index script
%
ave_hdr_offsets_v;
[fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid(1));

%
% Copy source header to destination file
%

if nfiles == 1
  ave_chdr = chdr;
  for c=1:fhdr(NCells)
    ave_chdr(c,NObs) = 1;
  end;
else
  ave_chdr = chdr;
end;
wt_ave_hdr_v(destfid,fhdr,ave_chdr,ename,cnames,fcom,ftext);

%
%  Determine the size of the input data array
%

nchan = fhdr(NChan);
npts = chdr(1,NPoints);
ncells = fhdr(NCells);

nsub = 0;

for ifile = 1:nfiles
  [fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, ...
	coff]=rd_egis_hdr_v(srcfid(ifile));
  nsub = nsub + chdr(1,NObs);
end

avgdata = zeros(npts,nchan,ncells,nsub); % Create empty output array

%
% Check to make sure we have bad chan info for each subject
%

if bad_chan & (size(bad_chan,1) ~= nsub)
  error(sprintf('Number of rows in bad chan list (%d) does not match number of subjects (%d)', size(bad_chan,1),nsub))
end

%
% Read the data from all files
%

curr_sub = 0;

for ifile=1:nfiles
  %
  % Read average file header
  %
  
  [fhdr,chdr,ename,czeros,cgains,cnames,fcom,ftext, coff]=rd_egis_hdr_v(srcfid(ifile));

  disp(sprintf('Loading data from %s', infname(ifile,:)))

  %
  % Begin looping through cells
  %

  curr_sub = curr_sub+1;

  for c=1:fhdr(NCells)
    
    %
    %  Loop through subject averages
    %

    nobs = chdr(c,NObs);
    good_obs = find(~ismember(1:nobs, bad_obs));
    ngood = length(good_obs);
    
    for t=1:ngood
  
      if t > 1
	curr_sub = curr_sub + 1;
      end
      
      %
      %  Read in the data from file
      %

      curr_obs = good_obs(t);
      avgdata(:,:,c,curr_sub) = rd_onetr_allch(srcfid(ifile), coff(c), curr_obs, ...
	  fhdr(NChan), chdr(c, NPoints));
      
      if ap.detrend
	disp(sprintf('  Detrending subject %d, cell %d ...', curr_sub, c))
	avgdata(:,:,c,curr_sub) = detrend(avgdata(:,:,c,curr_sub));
      end
      
    end 	%for t=1:chdr(c,NObs)
  end 	%for c=1:fhdr(NCells)
end	%for ifile =1:nfiles

disp(sprintf('Loaded data for %d subjects', nsub))

%
% Construct bad channel masks
%

mask = ones(npts,nchan,ncells,nsub);

for isub = 1:nsub
  
  % Add all channels with all zeros for any given cell and subject to mask
  for icell = 1:ncells
    mask(:,:,icell,isub) = repmat(any(avgdata(:,:,icell,isub)),npts,1);
%    disp(sprintf('Subject %d, cell %d badchans:', isub,icell))
%    disp(sprintf('%d ', find(~any(avgdata(:,:,icell,isub)))))
  end
  
  % Take into account bad channels passed as a parameter
  badidx = bad_chan(isub,find(bad_chan(isub,:)));
  mask(:,badidx,:,isub) = 0;
  
end  % for isub=

%
%  Compute z-score on a subject by subject basis if desired
%  Base the z-score on all good data from that subject
%

if ap.compute_z
  disp('Computing z-scores ...')
  for isub = 1:nsub
    disp(sprintf('  Subject %d', isub))
    tmpdata = avgdata(:,:,:,isub);
    tmpgood = find(mask(:,:,:,isub));
    mu = mean(tmpdata(tmpgood))
    sigma = std(tmpdata(tmpgood))
    tmpdata(tmpgood) = (tmpdata(tmpgood)-mu)/sigma;
    avgdata(:,:,:,isub) = tmpdata;
  end
  
  % Multiply by scaling factor
  avgdata= avgdata * 500;
  
end % if ap.compute_z

%
%  Compute the average
%

disp('Computing average')
gave_data = sum(avgdata,4) ./ sum(mask,4);

%
% Write data to .gave file
%

precision = 'integer*2';
disp(sprintf('Writing data to %s', outfname)) 
for c = 1:fhdr(NCells)
  fwrite(destfid, gave_data(:,:,c)', precision);
end;

disp('Finished running make_gave.');
fclose('all');

