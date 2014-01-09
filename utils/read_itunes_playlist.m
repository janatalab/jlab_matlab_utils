function data_st = read_itunes_playlist(fname)
% Reads a playlist text file generated by iTunes and returns an
% Ensemble-style data struct

% 25Dec2013 Petr Janata
logfid = 1;

if nargin < 1
  error('Filename containing playlist not specified');
elseif ~exist(fname,'file')
  error('File not found: %s', fname);
end

% Initialize the output data struct
data_st = ensemble_init_data_struct;
data_st.type = 'playlist';

data_st.meta.filename = fname;

% Open the file
fid = fopen(fname,'rt');

fprintf(logfid,'Reading playlist from: %s\n', fname);

% Read the headerline
hdrline = fgetl(fid);
vars = regexp(hdrline,'\t','split'); % Parse the variable names
vars = lower(regexprep(vars,'\s','_')); % Regularlize the varibale names
nvars = length(vars);

data_st.vars = vars; % copy variable names to output structure
cols = set_var_col_const(vars);

% Initialize the cell array
data_st.data = cell(1,nvars);

% Read in the rest of the playlist 
[data_st.data(:)] = textscan(fid, repmat('%s',1,nvars), 'delimiter','\t');

% Convert variables to desired types. Keep them as strings by default
vartypes.numeric = {'size','disc_number','disc_count','track_number','track_count', ...
  'bit_rate','sample_rate','plays','year','time','skips'};
vartypes.date = {'date_modified','date_added','last_played','last_skipped'};

for ivar = 1:nvars
  currVar = vars{ivar};
  
  % Get rows with missing values
  emptyMask = cellfun('isempty',data_st.data{ivar});
  switch currVar
    case vartypes.date
      tmpdatenum = nan(size(data_st.data{ivar}));
      
      % Replace missing values with bogus datenum
      [data_st.data{ivar}{emptyMask}] = deal(datestr(0));
      
      % Convert the date string to a datenum. Have to do this in two steps
      % because of datenum method
      if any(emptyMask)
        tmpdatenum(emptyMask) = datenum(data_st.data{ivar}(emptyMask));
      end
      if any(~emptyMask)
        tmpdatenum(~emptyMask) = datenum(data_st.data{ivar}(~emptyMask));
      end
      
      data_st.data{ivar} = tmpdatenum;
      
    case vartypes.numeric
      % Replace empty values with NaNs
      [data_st.data{ivar}{emptyMask}] = deal('NaN');
      
      % Convert to numeric
      data_st.data{ivar} = cellfun(@str2num, data_st.data{ivar});
      
      
  end % switch
end % for ivar=

% Close the file
fclose(fid);

return