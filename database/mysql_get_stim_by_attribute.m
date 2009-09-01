function [data, vars] = mysql_get_stim_by_attribute(varargin)
% Returns a list of stimuli associated with a given list of attributes
% [data, vars] = mysql_get_stim_by_attribute(varargin);
%
% Given a list of attribute names, this returns the list of stimuli that are associated
% with those attributes.
%
% Input arguments:
%   'name','attrib_name' - a string or cell array of strings containing the
%                          attribute names to be matched.
%   'conn_id' - database connection ID to use
%   'host' - database host
%   'database' - database on host to use
%   'params' - a structure containing the above parameters as fields
%   'extract_vars' | 'extract_flds' - specify column names to return, default=all
%       NOTE: all extract_vars must have the table specified (tbl.field,
%       not just 'field') ... due to the specialization of this function to
%       stim selection, only the tables specified in the 'tables' cell
%       array are valid.

% 07/12/07 Petr Janata
% 10/09/07 FB - general support to specify return variables, including support
% for when you pass '*'. if no extract_vars are passed, the default set will be used

data = {};
vars = {};
inputargs = varargin;

% valid search tables
tables = {'stimulus','attribute','stimulus_x_attribute'};

% Check to see if params is one of the input argumentes
str_idxs = find(cellfun(@isstr,inputargs));
params_idx = strmatch('params',inputargs(str_idxs));

if ~isempty(params_idx)
  params = inputargs{str_idxs(params_idx)+1};
  rm_idxs = str_idxs(params_idx):str_idxs(params_idx)+1;
  inputargs(rm_idxs) = [];
end

% Parse the input arguments
narg = length(inputargs);
for iarg = 1:2:narg
  switch inputargs{iarg}
    case 'params'
      params = inputargs{iarg+1};

    case 'conn_id'
      params.conn_id = inputargs{iarg+1};
      
    case 'host'
      params.host = inputargs{iarg+1};
    
    case 'database'
      params.database = inputargs{iarg+1};
      
    case {'name','attrib_name'}
      params.attrib_names = check_cell(inputargs{iarg+1});
    
    case {'extract_flds','extract_vars'}
      params.extract_vars = check_cell(inputargs{iarg+1});
      
  end % switch 
end % for iarg
 
try attrib_names = check_cell(params.attrib_names); catch attrib_names = {}; end
if isempty(attrib_names)
  sprintf('no attribute names specified\n')
  return
end

% Check for connection to database
try conn_id = params.mysql.conn_id;
catch
    try conn_id = params.conn_id;
    catch
        tmp_conn_id = 1;
        conn_id = 0;
    end
end

if mysql_check_conn(conn_id)
  try host = params.host; catch host = []; end
  try database = params.database; catch database = []; end
  mysql_make_conn(host,database,conn_id);
end

% Form the SQL query
if iscell(attrib_names)
  attrib_name_str = sprintf('"%s"', cell2str(attrib_names,'","'));
else
  attrib_name_str = sprintf('"%s"', attrib_names);
end

% if extract_vars are specified, use mysql_describe_table to make sure the
% requested variables exist as columns in the table. if no extract_vars are
% specified, pass * .. ?

if (~isfield(params,'extract_vars') || isempty(params.extract_vars))

    % default fields
    extract_vars = {'stimulus.stimulus_id','stimulus.name','attribute.name'};

% not sure if this is all necessary, or if I have covered all of the bases,
% to not have an error ...
elseif (isequal(params.extract_vars,{'*'}) || isequal(params.extract_vars,'*'))
    
    extract_vars = {};
    % get all fields from all tables, append table stub, populate extract_vars
    for itbls = 1:length(tables)
        % create mysql table stub
        tbl_stub = sprintf('%s.',tables{itbls});

        tbl_desc = mysql_describe_table(tables{itbls},conn_id);
        tbl_desc.flds = strcat(tbl_stub,tbl_desc.flds);
        extract_vars = [extract_vars transpose(tbl_desc.flds)];
    end
    
else

    % validate extract_vars members
    extract_vars = check_cell(params.extract_vars);
    tbl_members = {};
    valid = 1;
    for itbls = 1:length(tables)
        
        % create mysql table stub
        tbl_stub = sprintf('%s.',tables{itbls});

        % identify extract_vars indices for this table
        tbl_idxs = strmatch(tbl_stub,extract_vars);
        extract_vars{tbl_idxs};

        if (length(tbl_idxs) > 0)
            % track # of extract_vars members belonging to this table
            % get vars from table, validate
            tbl_desc = mysql_describe_table(tables{itbls},conn_id);
            for tidx = 1:length(tbl_idxs)
                tbl_members = [tbl_members extract_vars{tbl_idxs(tidx)}];
                candidate = strrep(extract_vars{tbl_idxs(tidx)},tbl_stub,'');
                if(~all(ismember(candidate,tbl_desc.flds)))
                    sprintf('%s does not exist within %s',candidate,tables{itbls})
                    valid = 0;
                end
            end
        end
    end

    % if we didn't account for all extract_vars members, then some of those
    % members did not map to one of the 'tables'
    if (length(tbl_members) ~= length(extract_vars))
        missing = setdiff(params.extract_vars,tbl_members);
        for midx = 1:length(missing)
            sprintf('%s is not within the allowed search tables\n',missing{midx})
        end
        valid = 0;
    end

    % if any invalid extract_vars were passed, do not continue
    if (~valid) return; end
        
end

extract_vars_str = cell2str(extract_vars,',');

mysql_str = sprintf(['SELECT %s FROM stimulus ' ...
      'LEFT JOIN stimulus_x_attribute USING (stimulus_id) ' ...
      'LEFT JOIN attribute USING (attribute_id) ' ...
      'WHERE attribute.name IN (%s);'], extract_vars_str, attrib_name_str);
									 
[data{1:length(extract_vars)}] = mysql(conn_id, mysql_str);

% If variables to extract contain periods, we need to remove those
for ivar = 1:length(extract_vars)
  varname =  extract_vars{ivar};
%   vars{ivar} = varname(findstr(varname,'.')+1:end);
%   varsrc = varname(1:findstr(varname,'.')-1); % for varsrc support,
%   varsrc must be added to the returns and caught upon return!!!!!!
  vars{ivar} = strrep(varname,'.','__');
end

if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
end % mysql_get_stim_by_attribute

function var = check_cell(var)
  if ~iscell(var)
    var = {var};
  end
end


