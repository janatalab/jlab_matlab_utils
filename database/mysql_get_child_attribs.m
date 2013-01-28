function [attribIDs, attribNames] = mysql_get_child_attribs(parentAttribs, params)
% Returns child attribute IDs and names associated with one or more parent
% attributes (specified either by name or ID) specified in the input. The
% attributes are found in the attribute_x_attribute table
%
% [attribIDs, attribNames] = mysql_get_child_attribs(parentAttribs, params);
%
% INPUTS: 
%    parentAttribs - can be a string, cell array of strings, or numeric
%                    vector specifying the name(s) or IDs, respectively of
%                    the parent attribute(s)
%    params - structure containing, most importantly, a mysql
%             field/structure with database information
%
% OUTPUTS:
%    attribIDs - a vector of attribute IDs
%    attribNames - a cell array of strings containing the child attribute
%                  names

% 27Jan2013 Petr Janata

attribIDs = [];
attribNames = {};

% Make sure our database connection is intact
if ~isfield(params,'mysql')
  error('params.mysql is required')
else
  conn_id = mysql_make_conn(params.mysql);
end

% Handle the input 
if ischar(parentAttribs)
  parentAttribs = {parentAttribs};
end

if isnumeric(parentAttribs)
  attribStr = sprintf('%d,',parentAttribs);
  attribStr(end) = '';
  mysql_str = sprintf('SELECT attribute_id FROM attribute WHERE attribute_id IN (%s) AND class="parent_category";', attribStr);
  parentAttribIDs = mysql(conn_id, mysql_str);
else
  attribStr = sprintf('"%s",', parentAttribs{:});
  attribStr(end) = '';
  mysql_str = sprintf('SELECT attribute_id FROM attribute WHERE name IN (%s) AND class="parent_category";', attribStr);
  parentAttribIDs = mysql(conn_id, mysql_str);
end

if isempty(parentAttribIDs)
  fprintf('Could not find any parent_category attribute IDs matching: %s\n', attribStr);
  return
end

parentStr = sprintf('%d,', parentAttribIDs);
parentStr(end) = '';

% Now consult the attribute_x_attribute table to retrieve the children
mysql_str = sprintf(['SELECT a.attribute_id, a.name FROM ' ...
  'attribute AS a, attribute_x_attribute AS axa ' ...
  'WHERE a.attribute_id = axa.attribute_id_child AND axa.attribute_id_parent IN (%s);'], parentStr);
[attribIDs, attribNames] = mysql(conn_id, mysql_str);

return