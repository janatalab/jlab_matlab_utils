function [Y,FS,NBITS] = mp3read(FILE)
%MP3READ Read MP3 (".mp3") sound file.
%    Y=MP3READ(FILE) reads a MP3 file specified by the string FILE,
%    returning the sampled data in Y. Amplitude values are in the range [-1,+1].
% 
%    [Y,FS,NBITS]=MP3READ(FILE) returns the sample rate (FS) in Hertz
%    and the number of bits per sample (NBITS) used to encode the
%    data in the file.
% 
%    Supports two channel encoded data, with up to 16 bits per sample.
% 
%    See also MP3WRITE, WAVWRITE, AUREAD, AUWRITE.

%%%%%% Location of the ".exe" Files
s = which('mp3read.m');
ww = findstr('mp3read.m',s);
location = s(1:ww-2);
mpg123 = location;;
mp3info = location;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Data extraction from File using "mp3info.exe"%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q = samprate, u = #frames, r = bitrate, v = mpeg version (1/2/2.5)
% C = Copyright, e = emph, E = CRC, L = layer, O = orig, o = mono, p = pad
cmd = [mp3info,'\mp3info', ' -r a -p  "%Q %u %r %v * %C %e %E %L %O %o %p" "', FILE,'"'];
w = mysystem(cmd);
% Break into numerical and ascii parts by finding the delimiter (*)
starpos = findstr(w,'*');
nums = str2num(w(1:(starpos - 2)));
strs = tokenize(w((starpos+2):end));
SR = nums(1);                           %Sampling Rate
nframes = nums(2);                      %Number of Frames
nchans = 2 - strcmp(strs{6}, 'mono');   %Number of Channels
layer = length(strs{4});                %MPEG Layer
bitrate = nums(3)*1000;                 %MPEG Bitrate
mpgv = nums(4);                         %MPEG Version
%%%%Temporary file%%%%%%
tmpfile = ['temp.wav'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data Decoding  using "mpg123.exe"%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd = [mpg123,'\mpg123', ' -w ', tmpfile, ' ', '"',FILE,'"'];
w = mysystem(cmd);
% Load the data and delete temporary file
[Y,FS,NBITS] = wavread(tmpfile);
delete(tmpfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = dos(cmd);
if s ~= 0 
  error(['unable to execute ',cmd]);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
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
    
