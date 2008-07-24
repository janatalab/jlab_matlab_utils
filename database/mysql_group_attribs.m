function mysql_group_attribs(parent_attrib_name,attrib_names,conn_id)
% Groups attributes together under a parent attribute in the attribute_x_attribute table
%
% mysql_group_attribs(parent_attrib_name,attrib_names,conn_id)
%

% 08/18/05 Petr Janata

% Do some input parameter checking. 
min_arg = 2;
max_arg = 3;

msg = nargchk(min_arg,max_arg,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

% Connect to host with a temporary connection if necessary
try conn_id(1);
catch   
  tmp_conn_id = 1;
  mysql_make_conn;
  conn_id = 0;
end

%
% Check to see if the parent attribute exists in the attribute table
%

mysql_str = sprintf(['SELECT attribute_id FROM attribute ' ...
      ' WHERE attribute.name="%s";'], parent_attrib_name);
parent_id = mysql(conn_id,mysql_str);

if isempty(parent_id)
  fprintf('Creating new parent attribute in attribute table ...\n');
  mysql_str = sprintf('INSERT INTO attribute (class,name) VALUES ("parent_category","%s");', parent_attrib_name);
  status = mysql(conn_id, mysql_str);
end

%
% Get the child attribute IDs
%

child_str = sprintf('"%s",',attrib_names{:});
child_str(end) = [];
mysql_str = sprintf(['SELECT attribute_id FROM attribute ' ...
      'WHERE name in (%s);'], child_str);
child_ids = mysql(conn_id, mysql_str);

%
% Now insert a bunch of child/parent relationships into the
% attribute_x_attribute table
%

value_str = sprintf('("%d","%d"),', [child_ids(:) ones(length(child_ids),1)*parent_id]');
value_str(end) = [];

mysql_str = sprintf(['INSERT INTO attribute_x_attribute' ...
      ' (attribute_id_child,attribute_id_parent) VALUES %s;'], value_str);
status = mysql(conn_id, mysql_str);

%
% Close the mysql connection if this was a temporary opening of the database
%

if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
