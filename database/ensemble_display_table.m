function ensemble_display_table(data_st,params)
% Mechanism for formatting display of data in a table.
%
% ensemble_display_table(data_st,params)
%
% The data to be displayed should be stored in a cell array in variable called
% data.  
% Additional variables are:
%   column_labels
%   column_formats

% 05/01/2007 Petr Janata
% 05/04/2007 PJ - fixed handling of cell arrays of strings that are passed
%                 in as a "single" value
% 07/24/2009 PJ - added optional replacement of NaNs

% Check integrity of data struct

% Set column indexing variable
cols = set_var_col_const(data_st.vars);

% Check display data to make sure we don't have any matrices within individual
% columns
ncols = size(data_st.data{cols.data},2);
for icol = 1:ncols
  if prod(size(data_st.data{cols.data}{icol})) ~= ...
	length(data_st.data{cols.data}{icol})
    msgstr = sprintf('%s: Display data column contains a matrix', mfilename);
    error(msgstr);
  else
    data_st.data{cols.data}{icol} = data_st.data{cols.data}{icol}(:);
  end
end

% Get a file identifier to write the table to. Defaults to standard out.
fid = ensemble_init_fid(params);

% Print some header information
try add_datestamp = params.add_datestamp; catch add_datestamp = false; end
if add_datestamp
  fprintf(fid,'\nTable generated on: %s\n\n', datestr(datenum(now),0));
end

% Print the column names
colnames = data_st.data{cols.column_labels};
if length(colnames) ~= ncols
  msgstr = sprintf(['%s: Mismatch between number of data columns (%d) and number of ' ...
	'column labels (%d)', ncols, length(colnames)]);
  error(msgstr)
end

for icol = 1:ncols
  fprintf(fid,'%s', colnames{icol});
  if icol < ncols, fprintf(fid,'\t');
  else fprintf(fid,'\n');
  end
end

% See whether we are going to replace NaNs
try replace_nan = params.replace_nan; catch replace_nan = false; end
if replace_nan
  try nan_str = params.nan_str; catch nan_str = '.'; end
end

% Print rows
nrows = length(data_st.data{cols.data}{1});
for irow = 1:nrows
  for icol = 1:ncols
    format_str = data_st.data{cols.column_formats}{icol};
    if ~isempty(findstr(format_str,'s'))
      if iscell(data_st.data{cols.data}{icol}{irow})
        data_val = cell2str(data_st.data{cols.data}{icol}{irow},',');
      else
        data_val = data_st.data{cols.data}{icol}{irow};
      end
    else
      data_val = data_st.data{cols.data}{icol}(irow);
      if isnan(data_val) && replace_nan
        data_val = nan_str;
        format_str = '%s';
      end
    end
    fprintf(fid, format_str, data_val);
    if icol < ncols, fprintf(fid,'\t');
    else fprintf(fid,'\n');
    end
  end
end

if fid > 1
  fclose(fid);
end
