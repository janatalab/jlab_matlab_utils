function sublist = make_exp_sublist(expname,varargin)
% Generates a list of unique subject ids for a given experiment.
% sublist = make_exp_sublist(expname,varargin);
%
% Generates a list of all unique subject IDs contained in the repsonse table of
% the experiment specified in expname.
%
%  'conn_id' - database connection ID to use during the query
% 
% Optional arguments:
%  'exclude_subs' - list of subjects to remove from the list

% 06/15/10 PJ - sanitized mysql_make_conn

conn_id = [];
exclude_subs = {};

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg

  switch varargin{iarg}
    case {'exclude_subs'}
      exclude_subs = check_cell(varargin{iarg+1});
         
    case 'conn_id'
      conn_id = varargin{iarg+1};
      
    otherwise
      fprintf('make_exp_sublist: Unknown parameter: %s\', varargin{iarg});
  end % switch
end % iarg

if ~exist('conn_id','var') || isempty(conn_id)
  error('%s: Do not have a valid connection ID', mfilename);
end

% Get the experiment info
expinfo = mysql_get_expinfo(expname,'','',conn_id);

sublist = expinfo.subs.ids;
  
% Remove subjects on exclusion list
if ~isempty(exclude_subs)
  [sublist] = setdiff(sublist,exclude_subs);
  sublist = sublist';
end

% Close the connection if necessary
if ~conn_id
  mysql(conn_id,'close')
end

end

function var = check_cell(var)
  if ~iscell(var)
    var = {var};
  end
end
