function [Y,FS,NBITS,OPTS] = mp3read(FILE,N,MONO,DOWNSAMP)
% [Y,FS,NBITS,OPTS] = mp3read(FILE,N,MONO,DOWNSAMP)   Read MP3 audio file
%       FILE is an MP3 audio file name.  Optional scalar N limits the number 
%       of frames to read, or, if a two-element vector, specifies 
%       the first and last frames to read.
%       Optional flag MONO forces result to be mono if nonzero;
%       Optional factor DOWNSAMP performs downsampling by 2 or 4.
%       Y returns the audio data, one channel per column; FS is the 
%       sampling rate.  NBITS is the bit depth (always 16). 
%       OPTS.fmt is a format info string.
% 2003-07-20 dpwe@ee.columbia.edu  This version calls mpg123.
% 2004-08-31 Fixed to read whole files correctly
% 2004-09-08 Uses mp3info to get info about mp3 files too
% 2004-09-18 Reports all mp3info fields in OPTS.fmt; handles MPG2LSF sizes
%            + added MONO, DOWNSAMP flags, changed default behavior.
% 2005-04-11 added -r m for report of median VBR in mp3info (thx Merrie Morris)
% 2006-01-09 Petr Janata - fixed argument list to mpg123 (see below).  Allowed
%            passing N = [];
% 2011-07-01 PJ - mpg123 would sometimes return insufficient data when only
%            a partial read was requested, so a bit extra is being requested now
% 2014-08-20 PJ - added dynamic mpg123, mp3info location based on
%                 architecture
% $Header: /homes/dpwe/matlab/columbiafns/RCS/mp3read.m,v 1.7 2005/04/11 20:55:16 dpwe Exp $

if nargin < 2
  N = 0;
else
  if isempty(N)
    N = 0;
  elseif length(N) == 1
    % Specified N was upper limit
    N = [1 N];
  end
end
if nargin < 3
  forcemono = 0;
else
  forcemono = (MONO ~= 0);
end
if nargin < 4
  downsamp = 1;
else
  downsamp = DOWNSAMP;
end
if downsamp ~= 1 && downsamp ~= 2 && downsamp ~= 4
  error('DOWNSAMP can only be 1, 2, or 4');
end

%%%%%%% Hard-coded behavior
%% What factor should we downsample by?  (1 = no downsampling, 2 = half, 4=qtr)
%downsamp = 4;
%% Do we want to force the data to be single channel? (1 = yes, 0= keep orig)
%forcemono = 1;

%%%%%% Location of the binaries
arch = computer('arch');
if strcmp(arch,'win32')
   fpath = fileparts(which(mfilename));
   mp3_toolbox_path = fullfile(fpath,'mp3_toolbox');
   if exist(mp3_toolbox_path)
     mpg123 = fullfile(mp3_toolbox_path,'mpg123');
     mp3info = fullfile(mp3_toolbox_path,'mp3info');
   else
     error('Cannot locate mpg123 and mp3info!')
   end
else
  mpg123 = '/usr/local/bin/mpg123';
  mp3info = '/usr/local/bin/mp3info';
end

%%%%%% Constants
NBITS=16;

%%%%%% Probe file to find format, size, etc. using "mp3info" utility
cmd = [mp3info, ' -r m -p "%Q %u %r %v * %C %e %E %L %O %o %p" "', FILE,'"'];
% Q = samprate, u = #frames, r = bitrate, v = mpeg version (1/2/2.5)
% C = Copyright, e = emph, E = CRC, L = layer, O = orig, o = mono, p = pad
w = mysystem(cmd);
if isempty(w) && exist(FILE,'file')
    trylim = 20;
    tries = 0;
    pauseval = 2;
    fprintf(1,'mp3info error, retrying %d times\n',trylim);
    while isempty(w)
      tries = tries+1;
      pause(pauseval);
      fprintf(1,'attempt %d after %0.2fs pause\n',tries,pauseval);
      w = mysystem(cmd);
      if tries == trylim
          error('%s found, but %s is not loading it',FILE,mp3info)
      end
    end % while isempty(W
end % if isempty(w) & exist(FILE,'file
    
% Break into numerical and ascii parts by finding the delimiter we put in
starpos = strfind(w,'*');
nums = str2num(w(1:(starpos - 2)));
strs = tokenize(w((starpos+2):end));

SR = nums(1);
nframes = nums(2);
nchans = 2 - strcmp(strs{6}, 'mono');
layer = length(strs{4});
bitrate = nums(3)*1000;
mpgv = nums(4);
% Figure samples per frame, after
% http://board.mp3-tech.org/view.php3?bn=agora_mp3techorg&key=1019510889
if layer == 1
  smpspfrm = 384;
elseif SR < 32000 && layer ==3
  smpspfrm = 576;
  if mpgv == 1
    error('SR < 32000 but mpeg version = 1');
  end
else
  smpspfrm = 1152;
end

OPTS.fmt.mpgBitrate = bitrate;
OPTS.fmt.mpgVersion = mpgv;
% fields from wavread's OPTS
OPTS.fmt.nAvgBytesPerSec = bitrate/8;
OPTS.fmt.nSamplesPerSec = SR;
OPTS.fmt.nChannels = nchans;
OPTS.fmt.nBlockAlign = smpspfrm/SR*bitrate/8;
OPTS.fmt.nBitsPerSample = NBITS;
OPTS.fmt.mpgNFrames = nframes;
OPTS.fmt.mpgCopyright = strs{1};
OPTS.fmt.mpgEmphasis = strs{2};
OPTS.fmt.mpgCRC = strs{3};
OPTS.fmt.mpgLayer = strs{4};
OPTS.fmt.mpgOriginal = strs{5};
OPTS.fmt.mpgChanmode = strs{6};
OPTS.fmt.mpgPad = strs{7};
OPTS.fmt.mpgSampsPerFrame = smpspfrm;

if SR == 16000 && downsamp == 4
  error('mpg123 will not downsample 16 kHz files by 4 (only 2)');
end

% 01/09/06 Petr Janata - changed downsampstr to be empty if downsamp rate is
% 1. This prevents only one channel being written when -1 is used to designate
% single channel writing in mpg123.

if downsamp == 1
  downsampstr = '';
else
  downsampstr = [' -',num2str(downsamp)];
end
FS = SR/downsamp;

if forcemono == 1
  nchans = 1;
  chansstr = ' -m';
else
  chansstr = '';
end

% Size-reading version
if strcmp(N,'size') == 1
   Y = [floor(smpspfrm*nframes/downsamp), nchans];
else

  % Temporary file to use
  tmpstub = sprintf('tmp%s.wav',num2str(round(1000*rand(1)))); 
  if strcmp(computer,'PCWIN')
    tmpdir = fullfile(FILE(1:regexp(FILE,filesep)-1),'TMP');
    if ~exist(tmpdir,'dir')
      mkdir(tmpdir)
    end
    tmpfile = fullfile(tmpdir,tmpstub);
  else
    tmpfile = fullfile('/tmp',tmpstub);
  end

  skipx = 0;
  skipblks = 0;
  skipstr = '';
  sttfrm = N(1)-1;
  if sttfrm > 0
    skipblks = floor(sttfrm*downsamp/smpspfrm);
    skipx = sttfrm - (skipblks*smpspfrm/downsamp);
    skipstr = [' -k ', num2str(skipblks)];
  end

  lenstr = '';
  endfrm = -1;
  if length(N) > 1
    endfrm = N(2);
    if endfrm > sttfrm
      decblk = ceil(endfrm*downsamp/smpspfrm) - skipblks + 1; % 01Jul2011 PJ - added 1 to ensure sufficient length read
      lenstr = [' -n ', num2str(decblk)];
    end
  end

  % Run the decode
  cmd=[mpg123, downsampstr, chansstr, skipstr, lenstr, ' -q -w ', tmpfile, ' ', '"',FILE,'"'];
  %disp(cmd);
  w = mysystem(cmd);

  % Load the data
  [Y,SR] = wavread(tmpfile);

  % Delete tmp file
  if strcmp(computer,'PCWIN')
    delete(tmpfile);
  else
    mysystem(['rm ', tmpfile]);
  end
  
  % Select the desired part
  if endfrm > sttfrm
    Y = Y(skipx+[1:(endfrm-sttfrm)],:);
  elseif skipx > 0
    Y = Y((skipx+1):end,:);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0 
  error(['unable to execute ',cmd]);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
% Debug
%disp([cmd,' -> ','*',w,'*']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = tokenize(s)
% Break space-separated string into cell array of strings
% 2004-09-18 dpwe@ee.columbia.edu
a = [];
p = 1;
n = 1;
l = length(s);
nss = findstr([s(p:end),' '],' ');
for ns = nss
  % Skip initial spaces
  if ns == p
    p = p+1;
  else
    if p <= l
      a{n} = s(p:(ns-1));
      n = n+1;
      p = ns+1;
    end
  end
end
    
