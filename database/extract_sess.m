function [table] = extract_sess(resp_tbl, sessids, conn_id)
% Extracts session information for a given session.
% [table] = extract_sess(resp_tlb, sessids, conn_id);
%
% Extracts all data associated with a specific session or list of sessions from
% an experiment response table.
%
% Also retrieves associated question information.
%
% resp_tbl -- name of table to extract response data from. 
% sessids - vector of session IDs to extract.  If sessids is empty, all
%           information is extracted.
%

% 06/22/06 PJ - adapted from extract_resp_tbl.m

global SQL_HOST DATABASE

if nargin < 1
  fprintf('extract_sess: Name of response table not specified');
  table = struct([]);
  return
end

try sessids(1);
catch  sessids = [];
end

% Connect to default host and database if no connection ID is specified
try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
end

%
% Determine which forms are associated with this experiment
%

nid = length(sessids)
for iid = 1:nid
  curr_sessid = sessids(iid);
  
  sql_str = sprintf(['SELECT DISTINCT form.form_id, form_name FROM form,%s ' ...
	'WHERE form.form_id=%s.form_id AND %s.session_id=%d'], resp_tbl, resp_tbl, resp_tbl, curr_sessid);
  [form_ids,form_names] = mysql(conn_id,sql_str);  % get the form ids

%  [form_names, idxs] = unique(form_names);
%  form_ids = form_ids(idxs);
  
  table{iid}.forms.names = form_names;
  table{iid}.forms.ids = form_ids;

  %
  % Loop through all of the forms and extract the data
  %

  nforms = length(form_names);
  fprintf('Found %d forms\n', nforms);
  table{iid}.forms.data = cell(nforms,1);
  for iform = 1:nforms
    curr_form = form_names{iform};
    
    fprintf('Processing form: %s\n', curr_form);
    switch curr_form
      case 'PANAS Survey2'
	table{iid}.forms.data(iform) = extract_panas(expt, subids, resp_tbl, conn_id);
      otherwise
	fprintf('Using generic handling for form: %s\n', curr_form);
	table{iid}.forms.data(iform) = extract_formdata(resp_tbl, form_ids(iform), ...
	    'conn_id', conn_id, ...
	    'sessids', curr_sessid);
    end % switch curr_form
  end % for  iform=
end % for iid

if ~conn_id
  mysql(conn_id, 'close');
end

return
