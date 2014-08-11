function data_st = ensemble_csv2datastruct(fname)
% Loads a CSV file into an Ensemble data structure

% 09Aug2014 Petr Janata

data_st = ensemble_init_data_struct;

if ~exist(fname,'file')
  error('%s: File %s does not exist', mfilename, fname)
end
[~,fstub] = fileparts(fname);
data_st.name = fstub;
data_st.type = 'csvfile';

% Oopen the file
fid = fopen(fname,'rt');

% Read the headerline
cl = fgetl(fid);

% Parse the line
vars = regexp(cl,',','split');
nvars = length(vars);

% Sanitize the variable names, replacing whitespace with underscores
vars = regexprep(vars,'\s+','_');
data_st.vars = vars;

% Read the rest of the file
data  = cell(1,nvars);
numRows = 0;

while ~feof(fid)
  cl = fgetl(fid); % read the line
  numRows = numRows+1;
  
  % Parse the line, taking care to preserve commas in quotes
  pat = '(".+")|([^,]*)';
  tokens = regexp(cl,pat,'match');
  
  % Replace double quotes
  tokens = regexprep(tokens,'"','');
  ntok = length(tokens);
  
  % Make sure the number of tokens is equal to the number of variables
  if ntok ~= nvars
    error('%s: Number of tokens (%d) does not match number of variables (%d)', mfilename, ntok, nvars)
  end
  
  % Determine whether variables are numeric or not
  if numRows == 1
    varIsNumeric = ~isnan(str2double(tokens));
  end
  
  for itok = 1:ntok
    % See if we should convert the value of each token to numeric
    % Place the token into the data slot
    if varIsNumeric(itok)
      data{itok}(numRows,1) = str2double(tokens{itok});
    else
      data{itok}{numRows,1} = tokens{itok};
    end
    
  end % for itok

end % while ~feof(fid)

% Close the file
fclose(fid);

% Make sure that the number of variables matches the number of columns
ncols = size(data,2);
if nvars ~= ncols
  error('%s: Number of columns (%d) does not match number of variables (%d)',mfilename, ncols, nvars);
end

data_st.data = data;
return