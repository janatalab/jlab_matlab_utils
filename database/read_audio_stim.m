function [sig,fs,nbits,opts] = read_audio_stim(stimulus_id,varargin)
% Reads in the information for an mp3 audio stim file from the stimulus table.
%
% [sig,fs,nbits,opts] = read_audio_stim(stimulus_id,N,mono,downsamp,mysql_conn_id)
%
% Accepts a stimulus_id, and optional parameters: N (the max. number of frames
% to read), mono (set to 1 to convert output to mono, 0 if no conversion), and
% downsamp, the downsampling factor
%
% Optionally, a mysql connection ID (mysql_conn_id) can be passed
% for communicating with the database. If this is omitted, the
% function will still be able to connect with the database, but
% this connection will result in connection information printing on
% the matlab console each time the function is called.
%
%
% The output variables are the same as those for wavread: sig (the signal 
% vector or matrix), fs (the sampling rate), nbits (number of bits per sample),
% and opts (a format info string).
%
% July 8, 2005 Original Version S.T.
% Sept. 15, 2005 Modified by S.T. use mp3read matlab function instead of sox 
% Dec. 7, 2005 Mod. by S.T. read stims from nfs space rather than copy to localhost using scp
% May 6, 2010 S.T. - added support for wav files. Since wavread doesn't support
%                    downsampling, this function simply reports an error if downsampling was set
%                    to a value other than one and a wav file was specified
  
ensemble_globals;

if(nargin < 5)
  mysql_conn_id = 7;
else
  mysql_conn_id = varargin{4};
end

if(nargin == 0)
  error('Need a stimulus_id from the stimulus table. Type ''help read_audio_stim'' for more information.');
end

if(nargin < 5)
  mysql_conn_id = 7;
end

if(mysql(mysql_conn_id,'status'))
  mysql_make_conn('','',mysql_conn_id);
end

sql_choose_stim = sprintf('select location from stimulus where stimulus_id = %d',stimulus_id);
stim_location = mysql(mysql_conn_id,sql_choose_stim);


stim_location = char(stim_location);

%determine the file extension and use the appropriate sound file reader for
%that extension
[p,fstub,ext] = fileparts(stim_location);

switch(ext) 
 case '.mp3'
  sfreader = @mp3read;
 case '.wav'
  sfreader = @wavread;
 otherwise
  error('Did not recognize sound file name extension');
end



%read the mp3 file
switch(nargin)
 case 1
  N = [];
  mono = 0;
  downsamp = 1;
 case 2
  N = varargin{1};
  mono = 0;
  downsamp = 1;
 case 3
  N = varargin{1};
  mono = varargin{2};
  downsamp = 1;
 case {4,5}
  N = varargin{1};
  mono = varargin{2};
  if(isempty(varargin{3}))
    downsamp = 1;
  else
    downsamp = varargin{3};
  end
  
end


if((downsamp ~= 1) && (strcmp(ext,'.wav')))
  error('This function doesn''t support downsampling of wav files. Wav files must be downsampled in the calling function.');
end

if(isempty(N) && (downsamp == 1))
  [sig,fs,nbits,opts] = sfreader(fullfile(stimulus_root,stim_location));
elseif(~isempty(N) && (downsamp == 1))
  [sig,fs,nbits,opts] = sfreader(fullfile(stimulus_root,stim_location),N);
elseif(isempty(N) && (downsamp ~= 1))
  %only mp3read supports this call (downsamp ~=1)
  [sig,fs,nbits,opts] = sfreader(fullfile(stimulus_root,stim_location),N,0,downsamp);
elseif(~isempty(N) && (downsamp ~= 1))
  %only mp3read supports this call (downsamp ~=1)
  [sig,fs,nbits,opts] = sfreader(fullfile(stimulus_root,stim_location),N,0,downsamp);
end

if(mono)
  sig = mean(sig,2);
end

