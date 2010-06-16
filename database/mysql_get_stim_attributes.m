function [stiminfo, stim_fields] = mysql_get_stim_attributes(stim_crit_fld, stim_crit_vals, conn_id)
% Returns info from the stimulus table for a given stimulus.
% [stiminfo, stim_fields] = mysql_get_stim_attributes(stim_crit_fld,stim_crit_vals);
%
% 'conn_id' - connection ID to database - required

% 11/24/06 PJ modified to handle list of stim ids
% 06/15/10 PJ msyql_make_conn handling

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

[tbl.flds,tbl.types,tbl.null, tbl.key,tbl.default,tbl.extra] = ...
    mysql(conn_id,'DESCRIBE stimulus');

stim_fields = tbl.flds;

if ~iscell(stim_crit_vals) && ~isstr(stim_crit_vals)
  stim_crit_str = sprintf('%d,', stim_crit_vals);
  stim_crit_str = stim_crit_str(1:end-1);
elseif iscell(stim_crit_vals)
  stim_crit_str = sprintf('"%s",', stim_crit_vals{:});
  stim_crit_str = stim_crit_str(1:end-1);
elseif isstr(stim_crit_vals)
  stim_crit_str = sprintf('"%s"', stim_crit_vals);
else
  fprintf('Do not know how to handle stim_crit_vals\n');
  return
end

% for compatability sake
if ~iscell(stim_crit_fld)
  stim_crit_fld = {stim_crit_fld};
end

sql_str = sprintf(['SELECT * FROM stimulus WHERE stimulus.%s ' ...
      'IN (%s);'], stim_crit_fld{1}, stim_crit_str);
[stiminfo{1:length(stim_fields)}] = mysql(conn_id,sql_str);

if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
