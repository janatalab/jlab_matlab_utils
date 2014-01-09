function attrib_id = mysql_stimfld2attrib(varargin)
% Associates stimuli with an attribute
%
% attrib_id = mysql_stimfld2attrib(varargin)
%
% Input arguments are tag/value pairs
%
% Tags:
% 'attribute' - name of attribute to be matched in attribute table
% 'stimtbl_fld' - stimulus table field to match values in
% 'stimtbl_vals' - values to search for in stimulus table field
% 'conn_id' - mysql connection ID

% 03/31/09 PJ
% 07/18/09 PJ - extracted portion of code into mysql_check_attribute
% 21Nov2013 PJ - added checking and conversion of string to cell during fldval handling

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case 'conn_id'
      conn_id = varargin{iarg+1};
    case {'attrib','attribute'}
      attrib_tag = varargin{iarg+1};
    case 'stimtbl_fld'
      fldname = varargin{iarg+1};
    case 'stimtbl_vals'
      fldvals = varargin{iarg+1};
    otherwise
      fprintf('%s: Unknown input argument: %s\n', mfilename, varargin{iarg});
  end
end

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

params.conn_id = conn_id;
params.create = true;
attrib_id = mysql_check_attribute(attrib_tag, params);

%
% Get the list of stimulus IDs that correspond to the current stimulus set
%
switch fldname
  case {'name','description','playlist','artist','album','genre','location'}
    if ~iscell(fldvals)
      fldvals = {fldvals};
    end
    stim_str = sprintf('"%s"',cell2str(fldvals,'","'));
  otherwise
    stim_str = sprintf('%d,', fldvals);
    stim_str(end) = '';
end
mysql_str = sprintf('SELECT stimulus_id FROM stimulus WHERE %s IN (%s);', fldname, stim_str);
stim_ids = mysql(conn_id,mysql_str);

%
% Get a list of stimulus IDs that are already associated with the attribute
%

mysql_str = sprintf(['SELECT stimulus_id FROM stimulus_x_attribute ' ...
    'WHERE attribute_id=%d'], attrib_id);
existing_stim_ids = mysql(conn_id, mysql_str);

% Determine which stim_ids to insert
insert_stim_ids = setdiff(stim_ids, existing_stim_ids);
insert_stim_ids = insert_stim_ids(:);

% Insert the new stim IDs
if ~isempty(insert_stim_ids)
  value_str = sprintf('("%d","%d"),', ...
      [insert_stim_ids ones(size(insert_stim_ids))*attrib_id]');
  value_str(end) = [];
  mysql_str = sprintf(['INSERT INTO stimulus_x_attribute (stimulus_id,' ...
	' attribute_id) VALUES %s;'], value_str);
  mysql(conn_id, mysql_str);
else
  fprintf('No stimulus IDs to be inserted\n');
end

return
