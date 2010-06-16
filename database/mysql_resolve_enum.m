function [strs,ids] = mysql_resolve_enum(dfid, conn_id)
% returns strings associated with different levels of an enum variable.
% Will accept a vector of dfids.  A cell array of cell string arrays
% will be returned.  For each dfid, there will be a cell array of strings
% corresponding to the different enum values.
%
% strs = mysql_resolve_enum(dfid);
%
% conn_id - connection to database - required

% 08/27/05 Petr Janata
% 11/10/05 PJ - modified output format
% 06/15/10 PJ - mysql_make_conn sanitization


% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

id_str = sprintf('%d,', dfid);
id_str(end) = [];

mysql_str = sprintf('SELECT enum_values, data_format_id FROM data_format WHERE data_format_id IN (%s);', id_str);
[tmp_strs, ids] = mysql(conn_id,mysql_str);

ndfid = length(dfid);
[strs{1:ndfid}] = deal({});

% Sort into original order
for iid = 1:ndfid
  enum_str = tmp_strs{ismember(ids,dfid(iid))};

  % Tokenize the string
  ie = 0;
  while ~isempty(enum_str)
    ie=ie+1;
    [strs{iid}{ie}, enum_str] = strtok(enum_str,',');
  end
end

return