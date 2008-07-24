function [outlist] = genre_cleaner(genre_mappings,conn_id)
% Renames variants of a genre name in the stimulus table
% of the database to a single genre name.
%

% 08/18/05 Petr Janata

% Connect to host with a temporary connection if necessary
try conn_id(1);
catch   
  conn_id = 0;
  mysql_make_conn;
end

%
% Get the master genre_map
%
genre_map = master_genre_map;
known_genres = cat(2,genre_map{:,2});

%
% Get a list of the genres currently present in the database
% 

db_genres = mysql_unique_fld_vals('stimulus','genre',conn_id);
num_db_genres = length(db_genres);

% First, loop through the list and make sure it has a mapping
unmapped = zeros(1,num_db_genres);

for ig = 1:num_db_genres
  if isempty(strmatch(db_genres{ig},known_genres,'exact'))
    unmapped(ig) = 1;
  end
end % for ig=

if any(unmapped)
  unmap_idx = find(unmapped);
  
  fprintf('Found %d unmapped genres: %s\n', length(unmap_idx),cell2str(db_genres(unmap_idx),','));
  fprintf(['Please incorporate them into master_genre_map.m and run this' ...
	' script again.\nMake sure you take into account whitespace\n']);
  outlist = db_genres(unmap_idx);
  return
end

% Now, for each genre variant that was encountered, check to see what its
% parent genre is. If these don't match, update the genre values in the
% stimulus table
ntarg_genres = size(genre_map,1);
for dbg = 1:num_db_genres
  for tg = 1:ntarg_genres
    % The next loop is necessary due to a bug in Matlab's strmatch function
    % var_idx = strmatch(db_genres{dbg}, genre_map{tg,2},'exact');
    nvar = length(genre_map{tg,2});
    var_idx = [];
    for ivar = 1:nvar
      if strcmp(db_genres{dbg},genre_map{tg,2}{ivar})
	var_idx = ivar;
	break
      end
    end
    
    if ~isempty(var_idx)  % are we dealing with the right parent genre
      % do we already match the parent genre
      parent_genre = genre_map{tg,1};
      curr_genre = genre_map{tg,2}{var_idx};
      if ~strcmp(curr_genre,parent_genre) % Need to update the database
	fprintf('Updating %s with %s\n', curr_genre, parent_genre);
	mysql_str = sprintf(['UPDATE stimulus ' ...
	      'SET stimulus.genre = "%s" ' ...
	      'WHERE genre = "%s";'], parent_genre, curr_genre);
	mysql(conn_id,mysql_str);
      end % if ~strcmp(genre_map{tg,2}{var_idx}, genre_map{tg,1})
    end % if ~isempty(var_idx)
  end % for tg=
end % for ig=

%
% Close the mysql connection if this was a temporary opening of the database
%

if ~conn_id
  mysql(conn_id,'close');
end
