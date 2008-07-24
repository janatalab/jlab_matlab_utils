% Lists stimuli in our mysql database
%
% stimlist.m
%

DEFAULT_HOST = 'atonal.ucdavis.edu';
DEFAULT_DATABASE = 'ensemble_main';

WRITE_TO_FILE = 0;

txtfile = 'songlist.txt';
min_duration = 30;  % minimum song sample duration (in s)
max_duration = 33;  % maximum song sample duration

try
  host(1);
catch
  host = DEFAULT_HOST;
end;

try 
  database(1);
catch
  database = DEFAULT_DATABASE;
end

% Connect to host
conn_id = mysql_make_conn(host);

% Switch to the experiment database
mysql('use',database);

[song,artist,album,year,totaltime,genre] = mysql('select name,artist,album,year,total_time,genre from stim_audio');

% Close the mysql connection
mysql('close');

% Determine song durations

[y,m,d,h,mi,s] = datevec(totaltime);

good_songs = find((s >= min_duration) & (s < max_duration));

if WRITE_TO_FILE
  % Write the song information to a text file
  fid = fopen(txtfile,'wt');

  fprintf(fid,'Song\tArtist\tAlbum\tYear\n');

  nsongs = length(good_songs);
  for isong = 1:nsongs
    song_idx = good_songs(isong);
    fprintf(fid,'%s\t%s\t%s\t%d\n', song{song_idx}, artist{song_idx}, album{song_idx}, year(song_idx));
  end

  fclose(fid);
end