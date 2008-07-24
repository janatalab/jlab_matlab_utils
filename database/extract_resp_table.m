function [table] = extract_resp_table(expt, resp_tbl, subids, conn_id)
% Extracts all data from an experiment response table.
% [table] = extract_resp_table(expt, table_name, subids, conn_id);
%
% Extracts all data from an experiment response table.
% Also retrieves associated question information.
%
% expt -- the name of the response table, where response_<expt> is the name of
% the table in the mysql database.
%
% table_name -- name of table to extract response data from. Overrides naming heuristic.
%

% 03/07/05 PJ
% 01/03/07 PJ - removed dependence on extract_generic()

global SQL_HOST DATABASE

try subids(1);
catch  subids = {};
end

try resp_tbl(1);
catch resp_tbl = '';
end

% Connect to default host and database if no connection ID is specified
try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
end

% Create a string identifying the response table
if isempty(resp_tbl)
  resp_tbl = sprintf('response_%s', lower(expt));
end

%
% Determine which forms are associated with this experiment
%

% This SQL query should be made more efficient
sql_str = sprintf('select form.form_id, form_name from form,%s where form.form_id=%s.form_id', resp_tbl, resp_tbl);
[form_ids,form_names] = mysql(conn_id,sql_str);  % get the form ids

[form_names, idxs] = unique(form_names);
form_ids = form_ids(idxs);

table.form_names = form_names;

%
% Loop through all of the forms and extract the data
%

nforms = length(form_names);
fprintf('Found %d forms\n', nforms);
table.form_data = cell(nforms,1);
for iform = 1:nforms
  curr_form = form_names{iform};
  
  fprintf('Processing form: %s\n', curr_form);
  switch curr_form
    case 'PANAS Survey2'
      table.form_data{iform} = extract_panas(expt, subids, resp_tbl, conn_id);
    otherwise
      fprintf('Using generic handling for form: %s\n', curr_form);
      try
	table.form_data{iform} = extract_generic(resp_tbl, ...
	    form_ids(iform),[],subids,conn_id);
      catch
	table.form_data(iform) = ...
	    extract_formdata(resp_tbl,form_ids(iform),'subids',subids,'conn_id',conn_id);
      end
  end % switch curr_form
end % for  iform=

if ~conn_id
  mysql(conn_id, 'close');
end

return
