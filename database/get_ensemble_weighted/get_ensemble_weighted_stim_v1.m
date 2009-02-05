function [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim(varargin)
% [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim(varargin);
%
% Returns the ID of a stimulus to present as a string.
%
% This is intended to be a very generic script that to manage stimulus
% selection needs of Ensemble experiments.  Experiment specific parameters
% should be stored in a paramters file, e.g. vdc_prescreen_params_v1.m, the
% name of which is passed into the function following a params_file tag.
%

script_version = 'get_ensemble_weighted_stim_v1.0';

% Parse input parameters
for iarg = 1:2:nargin
  curr_arg = varargin{iarg};
  
  switch curr_arg
    case {'params','params_file','param_file'}
      param_file = varargin{iarg+1};
    
  end % swich curr_arg
end % for iarg

if ~exist('param_file') || isempty(param_file)
  param_file = 'vdc_globals';
end

%
% Declare persistent variables
%
persistent is_initialized 
persistent conn_id
persistent params
persistent stimids stimMeta
persistent A W weight_cols

%
% Handle initialization
%
if isempty(is_initialized)  
  % Get parameters from a parameters file that is passed in as an input argument
  params_fh = str2func(param_file);
  params = params_fh();
  
  %Initialize random number generator
  rand('state', sum(100*clock));

  % Establish a connection to the database
  conn_id = params.mysql.conn_id;
  if mysql_check_conn(conn_id)
    mysql_make_conn(params.mysql.host, params.mysql.database, conn_id);
  end
  
  % Get a list of possible stimulus IDs (based on attribute/tag values in
  % stimulus_x_attribute table)
  [stimdata,stimvars] = mysql_get_stim_by_attribute('params',params);
  stim_cols = set_var_col_const(stimvars);
  
  % Get list of unique stimuli
  stim_ids = unique(stimdata{stim_cols.stimulus_id});
  nstim_total = length(stim_ids);
  
  % Get stimulus information
  stimMeta = mysql_extract_metadata('table','stimulus', ...
    'stimulus_id',stim_ids, ...
    'conn_id', conn_id);

  % Determine weighting parameters that we will be using and set these to their
  % initial values
  if ~isfield(params,'weights')
    fprintf(['No weighting specified. Using defaults of initial stimulus' ...
	  ' probability.\n']);
    params.weights.init_prob = 1;
    
    % Check to see if we are sampling with or without replacement and create a
    % weighting parameter if we are sampling without
    if isfield(params.selection,'replacement')
      params.weights.unplayed = ~params.selection.replacement;
    end
  end
  
  % Generate a list of weighting parameters and the weighting vector
  weight_names = fieldnames(params.weights);
  weight_cols = set_var_col_const(weight_names);
  num_weights = length(weight_names);
  W = zeros(1, num_weights);
  for iw = 1:num_weights
    W(iw) = params.weights.(weight_names{iw});
  end
  
  % Create an empty attribute matrix, A
  A = zeros(nstim_total,num_weights);
  
  % Initialize all stimuli to being unplayed
  A(:,weight_cols.unplayed) = ones(nstim_total,1);
  
  % Set stimulus weights
  switch params.selection.type
    case 'equal_probability'
      A(:,weight_cols.init_prob) = ones(nstim_total,1)/nstim_total;
     
    case 'weighted'
      
      % Analyze stimuli already presented in order to get their ratings (attribute
      % scores) if we will use these ratings as a basis for stim selection

      % The idea here is that this should be a function handle to the
      % experiment specific function that handles apriori weighting.  The
      % apriori weighting component could be a nested or sub-function in the
      % experiment weighting function, in the event that there is trial by
      % trial weighting function updating.

  end
  
  is_initialized = 1;
  if ~params.return_stim_during_init
    return
  end
end % if isempty(is_initialized)

% Update attributes if necessary. This should be an experiment specific
% function handle.  Default is no updating.

% Update variable weights if necessary. Experiment specific function
% handle. Default is no updating


% Pass the attribute matrix A, which is an M x N matrix where,
%    M = number of stimuli and 
%    N = number of attributes
% along with an attribute weighting vector, to select_from_weighted() which
% draws randomly from the distribution of stimulus probabilities.
  
item_idx = select_from_weighted(A,W);

% Set the output variables
stimloc = stimMeta(item_idx).location;
stim_id = stimMeta(item_idx).stimulus_id;
stim_id_str = num2str(stim_id);

% Update any columns that need to be updated based on this selection
A(item_idx,weight_cols.unplayed) = 0;

return
