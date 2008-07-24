function mp3write(Y,FS,NBITS,ENCODING,MP3FILE)
%MP3WRITE Write MP3 (".mp3") sound file.
%    MP3WRITE(Y,FS,NBITS,ENCODING,MP3FILE) writes data Y to a MP3
%    file specified by the file name MP3FILE, with a sample rate
%    of FS Hz and with NBITS number of bits. Stereo data should 
%    be specified as a matrix with two columns. 
%    ENCODING must be specified as an integer number from 1 to 5
% 
%   1 = Fixed bit rate 128kbs encoding.
%   2 = Fixed bit rate jstereo 128kbs encoding, high quality (recommended).
%   3 = Average bit rate 112kbs encoding.
%   4 = Fast encode, low quality.
%   5 = Variable bitrate.
%
%    See also MP3READ, WAVREAD, WAVWRITE.

s = which('mp3write.m');
ww = findstr('mp3write.m',s);
lame = s(1:ww-2);
wavwrite(Y,FS,NBITS,strcat(lame,'\temp.wav'));
tmpfile = strcat(lame,'\temp.wav');
MP3FILE = strcat(pwd,'\',MP3FILE);
ENCODING =  num2str(ENCODING);
switch ENCODING
    case {'1'}
        cmd = [lame,'\lame', ' --quiet', ' ', tmpfile, ' ',MP3FILE];
    case {'2'}
        cmd = [lame,'\lame', ' --quiet', ' -b 128 ', tmpfile, ' ',MP3FILE];
    case {'3'}
        cmd = [lame,'\lame', ' --quiet', ' --abr 112 ', tmpfile, ' ',MP3FILE];
    case {'4'}
        cmd = [lame,'\lame', ' --quiet', ' -f ', tmpfile, ' ',MP3FILE];
    case {'5'}
        cmd = [lame,'\lame', ' --quiet', ' -h ', ' -V ', tmpfile, ' ',MP3FILE];
    otherwise
        error('Encoding parameters not suported') 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data Encoding  using "Lame.exe"%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dos(cmd);
% Delete temporary file
delete(tmpfile);