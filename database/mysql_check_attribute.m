function attrib_id = mysql_check_attribute(attrib_tag, params)
% Checks to see whether an attribute exists in the attribute table and
% optionally creates it if necessary.
%
% attrib_id = mysql_check_attribute(attrib_tag, params);
% 
% attrib_tag is the attribute name to look for.  attrib_tag can be a cell array
% of strings in which case each string is checked as a tag, and a vector of
% attribute IDs is returned.
%
% If create is true and the attribute is not found, the attribute is created.
% If multiple tags are searched and the creation parameter is set to false, any
% missing attributes are indicated by NaNs in the attribute ID vector.

% 07/18/09 Petr Janata

% Check for connection to database
try 
  conn_id = params.conn_id;
  tmp_conn_id = 0;
catch
  tmp_conn_id = 1;
  conn_id = 0;
end

if mysql_check_conn(conn_id)
  try host = params.host; catch host = []; end
  try database = params.database; catch database = []; end
  mysql_make_conn(host,database,conn_id);
end

% Check to see if we are creating missing attributes
try create = params.create; catch create=false; end

% Check to see if attrib_tag is a cell array
if ~iscell(attrib_tag)
  attrib_tag = {attrib_tag};
end

num_attrib = length(attrib_tag);

for iattrib = 1:num_attrib

  % Check to see if the attribute tag already exists
  mysql_str = sprintf('SELECT attribute_id FROM attribute WHERE name="%s";', attrib_tag{iattrib});
  [curr_id] = mysql(conn_id, mysql_str);

  if isempty(curr_id) && create
    % Create the attribute if necessary
    mysql_str = sprintf(['INSERT INTO attribute (name,class) ' ...
	  'VALUES ("%s","stim_set");'], attrib_tag);
    mysql(conn_id,mysql_str);

    % Get the attribute ID
    mysql_str = sprintf('SELECT attribute_id FROM attribute WHERE name="%s";', attrib_tag);
    [attrib_id(iattrib)] = mysql(conn_id, mysql_str);
    fprintf('Created attribute (%s), ID=%d\n', attrib_tag, attrib_id(iattrib));
  elseif isempty(curr_id) && ~create
    if num_attrib == 1
      attrib_id = [];
    else
      attrib_id(iattrib) = NaN;
    end
  else
    attrib_id(iattrib) = curr_id;
    fprintf('Found attribute (%s), ID=%d\n', attrib_tag, attrib_id(iattrib));
  end
end

% Close connection if necessary
if tmp_conn_id
  mysql(conn_id,'close');
end
