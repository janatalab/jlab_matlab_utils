function [stiminfo, stim_fields] = mysql_get_expstims(resp_tbl,varargin)
% Retrieves auditory stimuli associated with a given response table.
% [stiminfo] = mysql_get_expstims(resp_tbl,varargin)
%
% Retrieves list of auditory stimuli (and associated info) of all stimuli that are
% associated with a response table.

% 12/29/07 Petr Janata

stiminfo = {};

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case 'conn_id'
      conn_id = varargin{iarg+1};
    otherwise
      fprintf('mysql_get_expstims: Unknown input argument: %s\n', varargin{iarg});
  end
end

% Check for connection to database
try conn_id(1);
catch   
  tmp_conn_id = 1;
  mysql_make_conn;
  conn_id = 0;
end

% Get list of stimulus IDs that are in the response table
mysql_str = sprintf(['SELECT DISTINCT stimulus_id FROM %s ' ...
      'WHERE stimulus_id IS NOT NULL;'], resp_tbl);
[stimids] = mysql(conn_id,mysql_str);

[stiminfo, stim_fields] = mysql_get_stim_attributes('stimulus_id',stimids,conn_id);

if (exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end

