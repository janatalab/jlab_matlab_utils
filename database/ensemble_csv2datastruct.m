function data_st = ensemble_csv2datastruct(in_st,params)
% Loads a CSV file into an Ensemble data structure
%
% USAGE:
% Either 
%    data_st = ensemble_csv2datastruct(fname)
% or
%    data_st = ensemble_csv2datastruct([],params)
%    where there is a field params.fname
%
% KNOWN ISSUES:
%   Fields containing strings with more than one comma will not parse
%   correctly. Need to improve the regexp statement.

% 09Aug2014 Petr Janata
% 26Jan2014 PJ - made compatible with ensemble_jobman
% 08Feb2015 PJ - parsing of input lines is now better able to handle empty
%                values, while retaining ability to preserve commas in
%                quotes

if nargin < 2
  if ischar(in_st)
    fname = in_st;
  else
    error('%s: argument must be string if only one argument is passed in', mfilename)
  end
  params = struct;
elseif nargin == 2
  if ~isfield(params, 'fname')
    error('%s: name of file to load must be provided in fname field in 2nd argument', mfilename)
  else
    fname = params.fname;
  end
end

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

if isfield(params,'USE_MATCH_HEURISTIC')
  USE_MATCH_HEURISTIC = params.USE_MATCH_HEURISTIC;
else
  USE_MATCH_HEURISTIC = 0;
end

while ~feof(fid)
  cl = fgetl(fid); % read the line
  numRows = numRows+1;
  
  
  % Parse the line, taking care to preserve commas in quotes
  % The problem with this approach is that if there is a missing value, the
  % number of tokens is less than the number of columns
  if USE_MATCH_HEURISTIC
    pat = '(".+")|([^,]*)'; %  '(".+")|([^,]*)'
    tokens = regexp(cl,pat,'match');
  else
    % Split the string on commas, but not those occuring in quotes
    % KNOWN ISSUE: Will give an erroneous result if the string within quotes
    % contains more than one comma. Need a better lookaround assertion in
    % the regexp pattern.
    pat = '(?!(?=,[\w\d\s]+")),';
    tokens = regexp(cl,pat,'split');

  end
  
  % Replace double quotes
  tokens = regexprep(tokens,'"','');
  ntok = length(tokens);
  
  % Make sure the number of tokens is equal to the number of variables
  if ntok ~= nvars
    % Check whether the data for the last variable is simply empty. This
    % would be true if the last character in the current line is a comma
    if ntok == (nvars-1) && strcmp(',',cl(end))
      tokens{end+1} = ' ';
      ntok = length(tokens);
    else
      error('%s: Number of tokens (%d) does not match number of variables (%d)', mfilename, ntok, nvars)
    end
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