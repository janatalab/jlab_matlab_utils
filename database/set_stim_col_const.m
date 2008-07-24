function SD = set_stim_col_const(var_list);
% Figures out which data columns correspond to different database fields in the stimulus table
%
% var_list is a cell string containing field names
% SD is a structure with constants (column ids) assigned for the following
% fields:
%
% STIM_ID -- stimulus ID
% NAME -- stimulus name
% DESC -- stimulus description
% ARTIST -- artist
% ALBUM -- album
% GENRE -- genre
% FORMAT -- file format
% SIZE -- file size
% DURATION -- duration of file HH:MM:SS
% YEAR -- year of recording
% BIT_RATE -- compression bit rate
% SAMPRATE -- sampling rate
% SAMPSIZE -- number of bits per sample?
% CHANS -- number of channels
% FNAME -- path to file

nvars = length(var_list);

FD = struct(...
    'STIM_ID',[], ...
    'NAME',[], ...
    'DESC',[], ...
    'ARTIST',[], ...
    'ALBUM',[], ...
    'GENRE',[], ...
    'FORMAT',[], ...
    'SIZE',[], ...
    'DURATION',[], ...
    'YEAR',[], ...
    'BIT_RATE',[], ...
    'SAMPRATE',[], ...
    'SAMPSIZE', [], ...
    'CHANS', [], ...
    'FNAME', [] ...
    );

for ivar = 1:nvars
  switch var_list{ivar}
    case 'stimulus_id'
      SD.STIM_ID = ivar;
    
    case 'name'
      SD.NAME = ivar;
      
    case 'description'
      SD.DESC = ivar;
      
    case 'artist'
      SD.ARTIST = ivar;
      
    case 'album'
      SD.ALBUM = ivar;

    case 'genre'
      SD.GENRE = ivar;
    
    case 'file_format'
      SD.FORMAT = ivar;
    
    case 'size'
      SD.SIZE = ivar;
    
    case 'duration'
      SD.DURATION = ivar;
    
    case 'year'
      SD.YEAR = ivar;
    
    case 'compression_bit_rate'
      SD.BIT_RATE = ivar;
    
    case 'sample_rate'
      SD.SAMPRATE = ivar;
    
    case 'sample_size'
      SD.SAMPSIZE = ivar;
    
    case 'channels'
      SD.CHANS = ivar;
    
    case 'location'
      SD.FNAME = ivar;
  
  end % switch
end % for ivar

return