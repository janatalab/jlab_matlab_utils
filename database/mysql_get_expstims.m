function [stiminfo, stim_fields] = mysql_get_expstims(resp_tbl,varargin)
% Retrieves auditory stimuli associated with a given response table.
% [stiminfo] = mysql_get_expstims(resp_tbl,varargin)
%
% Retrieves list of auditory stimuli (and associated info) of all stimuli that are
% associated with a response table.
%
% 'conn_id' <conn_id> is REQUIRED

% 12/29/07 Petr Janata
% 06/15/10 PJ - sanitized mysql_make_conn

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

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

% Get list of stimulus IDs that are in the response table
mysql_str = sprintf(['SELECT DISTINCT stimulus_id FROM %s ' ...
      'WHERE stimulus_id IS NOT NULL;'], resp_tbl);
[stimids] = mysql(conn_id,mysql_str);

[stiminfo, stim_fields] = mysql_get_stim_attributes('stimulus_id',stimids,conn_id);

return


