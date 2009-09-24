function primary_key = mysql_table_key(table,conn_id)
%returns the primary key name of a mysql table. If a composite key, returns a cell array of strings
% primary_key = mysql_table_key(table,conn_id)
%
% Author(s):
% Sept 23, 2009 - Stefan Tomic, first version
  
  
try
  conn_id;
catch
  conn_id= 0;
end

mysql_check_conn(conn_id,'open');

[tableDescription{1:6}] = mysql(conn_id,sprintf('describe %s',table));

%Key type is fourth column. There doesn't seem to be a general way to reference this
keyType = tableDescription{4};

primaryIdxs = strmatch('PRI',keyType);

if(length(primaryIdxs) == 1)
  primary_key = tableDescription{1}{primaryIdxs};
else
  primary_key = tableDescription{1}(primaryIdxs);
end
