function [data,vars] = stim_ids_from_resp_tbl(params)

% returns all unique stimulus ids from a given response table(s)
% 
%   [data,vars] = stim_ids_from_resp_tbl(params)
% 
% REQUIRED
%   params.stim_ids_from_resp_tbl.tables - struct array identifying
%       response tables from which to collect stimulus ids, as well as
%       optional filtering criteria to be applied to the data from each
%       table
%       .name - response table name
%       .filt - filtering criteria to be applied to the given table
%   params.mysql (params.ensemble) - database connection information
%       .host, .database, .conn_id
% 
% FIXME: uses the slightly unconventional output format of
% mysql_get_stim_by_attribute, with two returns: data and vars, as opposed
% to one return of an ensemble data struct containing data and vars
% 

% FB 2009.05.11
% 06/15/10 PJ - sanitized connection handling

% initialize vars
data{1} = nan(0);
vars = {'stimulus__stimulus_id'};

% initialize mysql struct, check mysql connection
if isfield(params,'mysql')
  m = params.mysql;
elseif isfield(params,'ensemble')
  m = params.ensemble;
else
  warning('no database connection info provided, assuming ensemble_main');
  m = struct();
end

if ~isfield(m,'conn_id') || mysql(m.conn_id,'status')
  if ~all(isfield(m,{'host','database','user','passwd'}))
    error('%s: No valid connection or insufficient information to establish one', mfilename);
  else
    m.conn_id = 6;
    mysql_make_conn(m);
    tmp_conn_id = 1;
  end
else
  tmp_conn_id = 0;
end

% get response table information
if ~isfield(params,'stim_ids_from_resp_tbl') || ...
        ~isfield(params.stim_ids_from_resp_tbl,'tables') || ...
        ~isstruct(params.stim_ids_from_resp_tbl.tables) || ...
        ~isfield(params.stim_ids_from_resp_tbl.tables,'name')
  error('required parameters not found\n');
else
  t = params.stim_ids_from_resp_tbl.tables;
end

% iterate over tables, extract unique stimulus_ids
for i = 1:length(t)

  % get data from the given response table
  tbl = t(i).name;
  tdata = ensemble_init_data_struct();
  tbl_desc = mysql_describe_table(tbl,m.conn_id);
  tdata.vars = transpose(tbl_desc.flds);
  tcol = set_var_col_const(tdata.vars);
  
  mysql_str = sprintf('SELECT %s FROM %s WHERE stimulus_id IS NOT NULL',...
      cell2str(tdata.vars,', '),tbl);

  [tdata.data{1:length(tdata.vars)}] = mysql(m.conn_id, mysql_str);
  
  % filter this data?
  if isfield(t(i),'filt') && isstruct(t(i).filt)
    tdata = ensemble_filter(tdata,t(i).filt);
  end
  
  % get the unique stimulus IDs
  b = unique(tdata.data{tcol.stimulus_id});
  data{1} = union(b,data{1});
end

if tmp_conn_id
  mysql(m.conn_id,'close');
  m.conn_id = [];
end
