function [trial_id] = two_afc_select_basic(varargin)
%  Selects trials for a Two-Alternative Forced Choice task
%
% [trial_id] = two_afc_select_basic(varargin);
%
% Copyright (c) 2006-13 The Regents of the University of California, Davis campus. All Rights Reserved.
%
% Author: Petr Janata
% Edits by JR
% Questions? contact jrector@ucdavis.edu

script_version = 'two_afc_select_basic_v1.0';

trial_id = '';

% Parse input parameters
for iarg = 1:2:nargin
  curr_arg = varargin{iarg};
  
  switch curr_arg
      case {'params','params_file','param_file'}
          param_file = varargin{iarg+1};    
  end % swich curr_arg
end % for iarg

% Make sure we have necessary variables
if ~exist('param_file') || isempty(param_file)
  fprintf('No parameter file specified\n');
  return
end

persistent is_initialized 
persistent master_trial_list

%
% Handle initialization
%
if isempty(is_initialized)  
  master_trial_list = [];
  
  %Initialize random number generator
  rand('state', sum(100*clock));
  
  % Get parameters from a parameters file that is passed in as an input argument
  params_fh = str2func(param_file);
  params = params_fh();

  % determine experiment name for logging
  if (isfield(params,'ensemble') && isfield(params.ensemble,'expname') && ~isempty(params.ensemble.expname))
    expname = params.ensemble.expname;
  else
    expname = 'unknown';
  end

  % Establish a connection to the database
  conn_id = params.mysql.conn_id;
  if mysql_check_conn(conn_id)
    mysql_make_conn(params.mysql);
  end

  %
  % Identify which trials are part of this experiment using the
  % trial_x_attribute table
  %
  
  % Retrieve the attribute ID
  mysql_str = sprintf(['SELECT attribute_id FROM attribute ' ...
	'WHERE name="%s";'], params.attrib_name);
  attrib_id = mysql(conn_id, mysql_str);
  
  % Get the list of trials
  mysql_str = sprintf(['SELECT trial_id FROM trial_x_attribute ' ...
	'WHERE attribute_id=%d'], attrib_id);
  trial_ids = mysql(conn_id, mysql_str);
  
  trial_id_str = sprintf('"%d",', trial_ids);
  trial_id_str(end) = '';
  
  % Load the trial information
  varlist = {'trial_id','data_format_id','correct_response_enum','stimulus_id1','stimulus_id2'};
  numvars = length(varlist);
  cols = set_var_col_const(varlist);
  
  mysql_str = sprintf(['SELECT %s FROM trial WHERE trial_id IN (%s);'], ...
      cell2str(varlist,','), trial_id_str);
  [td{1:numvars}] = mysql(conn_id, mysql_str);
  
  % Check to see if we are dealing with homogeneous trial types
  dfid = unique(td{cols.data_format_id});
  if length(dfid) > 1
    fprintf('Too many data format IDs associated with trials\n');
    return
  end  
  
  % Add all available trials to the list if param is set accordingly
  if params.include_all_trials
    master_trial_list = td{cols.trial_id};
  else
    % implement logic here to determine which trials to include
  end

  % Shuffle the master trial list
  num_trials = length(master_trial_list);
  master_trial_list = master_trial_list(randperm(num_trials));

  fprintf('%s for experiment ''%s'': Initialized %d trials\n', mfilename, expname, num_trials);
  is_initialized = true;
end % if isempty(is_initialized)  

% Select the first trial ID from the list and eliminate it from the list
if ~isempty(master_trial_list)
  trial_id = master_trial_list(1);
  master_trial_list(1) = [];
  
  trial_id = num2str(trial_id);
else
  trial_id = 'end';
end
