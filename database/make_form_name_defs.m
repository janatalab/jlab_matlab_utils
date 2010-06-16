function [form_id_const, form_name_id_const_map, form_name_list] = make_form_name_defs(params)
% Generates the form_name_defs file if params.write2file is true. Otherwise it
% returns a set of variables that allow one to get form IDs by name.
%
% [form_id_const, form_name_id_const_map, form_name_list] = make_form_name_defs(params);
%
%  Constants are now fields of a form_name_const structure.  form_name_list is a
%  corresponding cell array that contains the original (Ensemble) name of the form
%  and the new constant name.  In order to get at the form ID field from the
%  original name, one can write form_name_const.(form_name_list{form_idx,2})

% 10/28/06 Petr Janata - updated from an earlier version.  
% 01/03/07 PJ - fixed to implement stated functionality and deal with bad
%          characters
% 10/15/07 PJ - Turned into a function, fixed handling of multiple databases
%          and added returning of the variables
% 12/11/08 PJ - Turned off writing of form_name_defs.m file by default
% 8/17/09 Stefan Tomic - added flag to only close db connection if
%                        the ID wasn't passed into the function
%                        (otherwise, this is a conn that we may want to use later)
% 06/15/10 PJ - sanitized mysql_make_conn

close_conn = 0;

if nargin && isfield(params,'ensemble')
  try CONN_ID = params.ensemble.conn_id; catch CONN_ID = 0; close_conn = 1; end
  try DATABASE_SCRIPT_PATH = params.paths.project_root; 
  catch DATABASE_SCRIPT_PATH = '.'; 
  end
  
  try write2file = params.write2file; catch write2file = 0; end
  
else
  CONN_ID = 0;
  DATABASE_SCRIPT_PATH = './';
  write2file = 0;
  close_conn = 1;
end

mysql_make_conn(params.ensemble);

mysql_str = sprintf('SELECT form_name, form_id FROM form');
[names, ids] = mysql(CONN_ID, mysql_str);

if write2file
  fname = fullfile(DATABASE_SCRIPT_PATH, 'form_name_defs.m');

  % Make a backup copy of the defs file if it already exists
  if exist(fname)
    unix_str = sprintf('mv form_name_defs.m form_name_defs.m.%s', datestr(datenum(now),29));
    unix(unix_str);
  end

  fid = fopen(fname, 'wt');
  if fid == -1
    error(sprintf('Could not open file: %s\n', fname))
  end

  fprintf('Writing %s\n', fname);

  fprintf(fid, ...
      ['%% Automatically generated mapping of form names to form numbers.\n%%\n%%'...
	' form_name_defs.m\n%%\n%% Generated on %s by make_form_name_defs.m\n\n'], datestr(datenum(now)));
end

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
  if write2file
    fprintf(fid,'form_id_const.%s = %d;\n', new_name, ids(iform));
  end
  form_id_const.(new_name) = ids(iform);
end

% Write a list of form name mappings
if write2file
  fprintf(fid, '\n\nform_name_id_const_map = { ...\n');
end

for iform = 1:nforms
  if write2file
    fprintf(fid,'''%s'', ''%s''; ...\n', strrep(form_name_list{iform,1},'''',''''''), form_name_list{iform,2});
  end
  form_name_id_const_map{iform,1} = ...
      strrep(form_name_list{iform,1},'''','''''');
  form_name_id_const_map{iform,2} = form_name_list{iform,2};
end
if write2file
  fprintf(fid, '};\n');
end

% Write the list of form names to a variable at the end of the file
if write2file
  fprintf(fid, '\n\nform_name_list = { ...\n');
end
for iform = 1:nforms
  if write2file
    fprintf(fid,'''%s''; ...\n', strrep(form_name_list{iform,1},'''',''''''));
  end
  new_form_name_list{iform,1} = strrep(form_name_list{iform,1},'''','''''');
end

if write2file
  fprintf(fid, '};\n');
  fclose(fid);
end

form_name_list = new_form_name_list;


% Close the mysql connection
if(close_conn)
  mysql('close');
end
