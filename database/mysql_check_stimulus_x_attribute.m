function out_st = mysql_check_stimulus_x_attribute(attrib_dict, params)
% Checks to see whether an attribute exists in the attribute table and
% optionally creates it if necessary.
%
% out_st = mysql_check_attribute(attrib_dict, params);
% 
% INPUT:
% 
% attrib_dict is a cell array of structures containing the following fields to match/insert:
%   'stimulus_id' - the ID of the stimulus in the stimulus table
%   'attribute_id' - the ID of the attribute in the attribute table
%   'attribute_value_double' - a numeric value that is optionally associated with the record
%   'attribute_value_text', - a string (maxlen=128) that is optionally associated with the record
%
% If attrib_dict is an array, a vector of attribute IDs is returned.
%
% params.conn_id - specify the mysql connection to use
% params.create - indicate whether attribute should be created if not
%                  present [default: false]
%
% If create is true and the attribute is not found, the attribute is created.
%
% OUTPUT:
% An Ensemble data struct

% 23May2015 Petr Janata - adapted from mysql_check_attribute()

% Check for connection to database
try 
  conn_id = params.conn_id;
  tmp_conn_id = 0;
catch
  tmp_conn_id = 1;
  conn_id = 0;
  params.conn_id = conn_id;
end

if mysql_check_conn(conn_id)
  mysql_make_conn(params);
end

out_st = ensemble_init_data_struct;

% Check to see if we are creating missing attributes
if isfield(params,'create')
  create = params.create; 
else
  create=false; 
end

if isstruct(attrib_dict) && length(attrib_dict)==1
  attrib_dict = {attrib_dict};
end

num_entries = length(attrib_dict);
tbl = mysql_describe_table('stimulus_x_attribute',conn_id);
vars = tbl.flds(:)';
data = cell(1,length(vars));
out_st.vars = vars;
out_st.data = data;

for ientry = 1:num_entries
  curr_st = ensemble_init_data_struct;
  curr_st.vars = vars;
  
  curr_entry = attrib_dict{ientry};
  flds = fieldnames(curr_entry);
  nflds = length(flds);
  
  % Construct the where clause for the query based on the information we
  % have
  insertColsStrCells = cell(1,nflds);
  whereStrCells = cell(1,nflds);
  valuesStrCells = cell(1,nflds);
  
  for ifld = 1:nflds
    currFld = flds{ifld};
    insertColsStrCells{ifld} = currFld;
    
    switch currFld
      case 'attribute_value_text'
        typstr = '%s';
      otherwise
        typstr = '%d';
    end
    fmtstr = sprintf('%%s="%s"', typstr);
    
    currVal = curr_entry.(flds{ifld});
    whereStrCells{ifld} = sprintf(fmtstr, flds{ifld}, currVal);
    valuesStrCells{ifld} = sprintf(['"' typstr '"'], currVal);
  end
  whereStr = cell2str(whereStrCells,' AND ');
  
  % Check to see if the attribute tag already exists
  mysql_str = sprintf('SELECT * FROM stimulus_x_attribute WHERE %s;', whereStr);
  [data{:}] = mysql(conn_id, mysql_str);
  foundEntries = ~isempty(data{1});

  if ~foundEntries && create
    % Create the attribute if necessary
    mysql_str = sprintf(['INSERT INTO stimulus_x_attribute (%s) ' ...
	  'VALUES (%s);'], cell2str(insertColsStrCells,','), cell2str(valuesStrCells,','));
    mysql(conn_id,mysql_str);

    % Get the attribute ID
    mysql_str = sprintf('SELECT * FROM stimulus_x_attribute WHERE %s;', whereStr);
    [data{:}] = mysql(conn_id, mysql_str);
    fprintf('Created entry: %s\n', whereStr);
  elseif ~foundEntries && ~create
    data = {};
  else
    fprintf('Found entry: %s\n', whereStr);
  end
  curr_st.data = data;
  
  % Append the datastruct
  out_st = ensemble_concat_datastruct({out_st, curr_st});
end

% Close connection if necessary
if tmp_conn_id
  mysql(conn_id,'close');
  params.conn_id = [];
end

return