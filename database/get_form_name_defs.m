function [form_id_const, form_name_id_const_map, form_name_list] = get_form_name_defs(params)
% Loads form name defs
% [form_id_const, form_name_id_const_map, form_name_list] = get_form_name_defs(params);
%
%  Constants are now fields of a form_name_const structure.  form_name_list is a
%  corresponding cell array that contains the original (Ensemble) name of the form
%  and the new constant name.  In order to get at the form ID field from the
%  original name, one can write form_name_const.(form_name_list{form_idx,2}).
%
%  This script is very similar to make_form_name_defs() with the exception that
%  it doesn't write an m-file (form_name_defs.m)

% 10/21/07 PJ - adapted from make_form_name_defs.m  The difference is that this
%               script doesn't write out a file.

if nargin && isfield(params,'ensemble')
  try HOST = params.ensemble.host; catch HOST = ''; end
  try DB = params.ensemble.database; catch DB = ''; end
  try CONN_ID = params.ensemble.conn_id; catch CONN_ID = 0; end
else
  HOST = '';
  DB = '';
  CONN_ID = 0;
end

mysql_make_conn(HOST,DB,CONN_ID);

mysql_str = sprintf('SELECT form_name, form_id FROM form');
[names, ids] = mysql(CONN_ID, mysql_str);

nforms = length(ids);
form_name_list = cell(nforms,2);
for iform = 1:nforms
  orig_name = names{iform};
  
  % Replace any characters in the name that won't work as variable names
  
  curr_name = regexprep(orig_name, {'-',' ','/',','}, '_');
  curr_name = regexprep(curr_name, {'''','(',')','\.'}, '');
  new_name = upper(curr_name);
  
  form_name_list{iform,1} = orig_name;
  form_name_list{iform,2} = new_name;
  form_id_const.(new_name) = ids(iform);
end

% Write a list of form name mappings
for iform = 1:nforms
  form_name_id_const_map{iform,1} = ...
      strrep(form_name_list{iform,1},'''','''''');
  form_name_id_const_map{iform,2} = form_name_list{iform,2};
end

% Write the list of form names to a variable at the end of the file
for iform = 1:nforms
  new_form_name_list{iform,1} = strrep(form_name_list{iform,1},'''','''''');
end

form_name_list = new_form_name_list;

% Close the mysql connection
mysql('close');
