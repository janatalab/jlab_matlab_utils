function [data,vars] = mysql_extract_data(varargin)
% This is a very generic routine for extracting data from a table.
%
% [data,vars] = mysql_extract_data(varargin);
%
% This is a very generic routine for extracting data from a table.
% It accepts an variable length argument list of tag/value pairs.
%
% 'table' - specifies name of table to extract data from
% {'extract_flds','extract_vars'} - a cell array of strings containing a list of fields to
%                  extract from the table. If this is not specified, all fields
%                  are extracted.
% {'order_by','sort_by'} - an optional field that specifies which fields should be used to
%              order the data
%
% 'conn_id' - the connection ID to the MySQL database to use. If this is left
%             empty, a temporary ID is created.
%
% All other tag/value pairs specify fields to be used for searching and the
% values to search for. With the exception of the following variables, the tag
% must match the field name in the table exactly.
%
% {'session_id','session_ids','sessid','sessids'} - map to 'session_id'
% {'subject_id','subject_ids','subid','subids'} - map to 'subject_id'
% {'form_id','form_ids','formid','formids'} - map to 'form_id'

% 01/04/07 Petr Janata
% 01/10/07 S.T. Added clause for handling empty search criteria
%               (i.e. no extract_flds or values are given)
% 03/12/08 PJ field names now enclosed in quotes. MySQL 5 compatibility


% Initialize some variables
fld.crit_flds = {};
fld.extract_flds = {};
fld.order_by = {};
crit_vals = {};
table = '';

data = {};
vars = {};

check_fld_vars = fieldnames(fld);

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg

  switch varargin{iarg}
    case {'experiment_id','experiment_ids','expid','expids'}
      fld.crit_flds{end+1} = 'experiment_id';
      crit_vals{end+1} = check_cell(varargin{iarg+1});
   
    case {'session_id','session_ids','sessid','sessids'}
      fld.crit_flds{end+1} = 'session_id';
      crit_vals{end+1} = check_cell(varargin{iarg+1});
   
    case {'subject_id','subject_ids','subid','subids'}
      fld.crit_flds{end+1} = 'subject_id';
      crit_vals{end+1} = check_cell(varargin{iarg+1});
   
    case {'form_id','form_ids','formid','formids'}
      fld.crit_flds{end+1} = 'form_id';
      crit_vals{end+1} = check_cell(varargin{iarg+1});
      
    case {'extract_flds','extract_vars'}
      fld.extract_flds = check_cell(varargin{iarg+1});
      
    case {'order_by','sort_by'}
      fld.order_by = check_cell(varargin{iarg+1});
      
    case {'table','resp_table'}
      table = varargin{iarg+1};
      
    case 'conn_id'
      conn_id = varargin{iarg+1};
      
    otherwise % assume a criterion field/value pair
      fld.crit_flds{end+1} =  varargin{iarg};
      crit_vals{end+1} = check_cell(varargin{iarg+1});
      
  end % switch
end % iarg

% Make sure a table was specified
if isempty(table)
  fprintf('mysql_extract_data: No table specified\n');
  return
end

% Check for connection to database
try conn_id(1);
catch 
  tmp_conn_id = 1;
  CONN_ID = 0;
  conn_id = mysql_make_conn([],[],CONN_ID);
end

% Make sure that all the fields we want to extract or sort by actually exist in the table
% we want to extract them from
tbl = mysql_describe_table(table,conn_id);

nfld_vars = length(check_fld_vars);
for ivar = 1:nfld_vars
  curr_fld = check_fld_vars{ivar};
  if strcmp(curr_fld,'extract_flds') && isempty(fld.extract_flds)
    fld.extract_flds = tbl.flds;
  end

  exist_mask = ismember(fld.(curr_fld),tbl.flds);
  bad_flds = find(~exist_mask);
  if ~isempty(bad_flds)
    for ibad = 1:length(bad_flds)
      fprintf('Field %s does not exist in table %s\n', fld.(curr_fld){bad_flds(ibad)}, table);
    end
    fld.(curr_fld)(bad_flds) = [];  % remove bad fields
  end % if ~isempty(bad_flds)
end % for ivar

% Prepare elements of the query
extract_vars_str = cell2str(fld.extract_flds,',');
extract_vars_str = sprintf('`%s`,',fld.extract_flds{:});
extract_vars_str(end) = [];

% Make the criterion string
ncrit = length(fld.crit_flds);
if(ncrit == 0)
  crit_str = [];
else
  crit_str = 'WHERE';
  for icrit = 1:ncrit
    if icrit > 1
      crit_str = [crit_str ' AND'];
    end
  
    if isstr(crit_vals{icrit}{1})
      crit_val_str = sprintf('"%s",', crit_vals{icrit}{:});
    else
      crit_val_str = sprintf('%d,', crit_vals{icrit}{:});
    end
    crit_val_str(end) = [];
    
    curr_crit_str = sprintf(' %s IN (%s)', fld.crit_flds{icrit},crit_val_str);
    
    crit_str = [crit_str curr_crit_str];
  end
end %else
  
  
sort_str = '';
if ~isempty(fld.order_by)
  sort_str = sprintf('ORDER BY %s', cell2str(fld.order_by,','));
end

% Construct the query 
mysql_str = sprintf(['SELECT %s FROM %s ' ...
      '%s %s'], extract_vars_str, table, crit_str, sort_str);

% Extract the data
[data{1:length(fld.extract_flds)}] = mysql(conn_id, mysql_str);

vars = fld.extract_flds;

if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
end % mysql_extract_data

function var = check_cell(var)
  if ~iscell(var)
    var = {var};
  end
end
