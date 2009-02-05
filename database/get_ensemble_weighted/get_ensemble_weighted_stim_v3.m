function [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim_v3(varargin)

% [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim_v3(varargin);
%
% Returns the ID of a stimulus to present as a string.
%
% This is intended to be a very generic script that to manage stimulus
% selection needs of Ensemble experiments.  Experiment specific parameters
% should be stored in a paramters file, e.g. vdc_prescreen_params_v1.m, the
% name of which is passed into the function following a params_file tag.
%
% FB Sep 27 2007 - adding functionality to weight stims, as well as to 
% dynamically update stim weighting ... this has been moved out of the
% 'is_initialized' conditional, so that it can be continuously updated
% 


script_version = 'get_ensemble_weighted_stim_v3.0';

% Parse input parameters
for iarg = 1:2:nargin
  curr_arg = varargin{iarg};
  
  switch curr_arg
      case {'params','params_file','param_file'}
          param_file = varargin{iarg+1};
      case {'subid','subject_id','subjectid'}
          subid = varargin{iarg+1};
    
  end % swich curr_arg
end % for iarg

if ~exist('param_file') || isempty(param_file)
  param_file = 'vdc_globals';
end

if ~exist('subid')
    sprintf('No Subid specified')
    stim_id_str = '';
    stim_id = '';
    stimloc = '';
    return
end

%
% Declare persistent variables
%
persistent is_initialized 
persistent conn_id
persistent params
persistent stimids stimMeta smCols
persistent A W weight_cols Wi prev_stim

%
% Handle initialization
%
if isempty(is_initialized)  
    
  % Get parameters from a parameters file that is passed in as an input argument
  params_fh = str2func(param_file);
  params = params_fh();
  params.subid = subid;
  %Initialize random number generator
  rand('state', sum(100*clock));

  % Establish a connection to the database
  conn_id = params.mysql.conn_id;
  if mysql_check_conn(conn_id)
    mysql_make_conn(params.mysql.host, params.mysql.database, conn_id);
  end
  
  % % % % 
  % identify this experiment, mark it for this subject
  
  % Get the ID of this experiment
  mysql_str = sprintf('SELECT experiment_id FROM experiment WHERE response_table="%s";',...
      params.resptbl_name);
  curr_expid = mysql(conn_id, mysql_str);

  %   mysql(conn_id,sprintf('use %s',params.mysql.database));
  % Insert the script name into the table
  insert_str = sprintf(['UPDATE %s SET %s.misc_info = ''param_file:%s'' ' ...
      'WHERE subject_id=''%s'' AND experiment_id=%d AND form_id=%d;'], ...
    params.resptbl_name,params.resptbl_name,params.version,...
    subid,curr_expid,params.targ_form_id);
  
  mysql(conn_id,insert_str);

  % if stim_compiler is defined, then run that function, else get the
  % attributes as they are defined within attrib_names, through
  % mysql_get_stim_by_attribute
  if (isfield(params,'fun') & isfield(params.fun,'stim_compiler') ...
          & isa(params.fun.stim_compiler,'function_handle'))
      [stimdata,stimvars] = params.fun.stim_compiler(params);
  else
      % Get a list of possible stimulus IDs (based on attribute/tag values in
      % stimulus_x_attribute table)
      [stimdata,stimvars] = mysql_get_stim_by_attribute('params',params);
  end
  
  stim_cols = set_var_col_const(stimvars);

  % Get list of unique stimuli
  stim_ids = unique(stimdata{stim_cols.stimulus__stimulus_id});
  nstim_total = length(stim_ids);
   
  % Get stimulus information
  [stimMeta.data,stimMeta.vars] = mysql_extract_data('table','stimulus', ...
    'stimulus_id',stim_ids, ...
    'conn_id', conn_id);

  smCols = set_var_col_const(stimMeta.vars);

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
  
  if (isfield(params.selection,'replacement') && ...
          ~params.selection.replacement)

      % Get info on this subject including the experiments they've been in
      sinfo = mysql_get_sinfo(subid,conn_id);
      sinfo.subid = subid;

      % % % %   
      % get stim this participant has already heard

      % Figure out if the subject was in previous versions of this
      % experiment, and get a list of stimuli that were presented to the subject
      prev_exp_idxs = find(ismember(sinfo.exp_names, params.prev_expts));
      prev_stim_ids = find_presented_stims(sinfo, prev_exp_idxs, conn_id);
      prev_stim_idxs = find(ismember(stimMeta.data{smCols.stimulus_id},...
          prev_stim_ids));
      A(prev_stim_idxs,:) = 0;
%       A(:,weight_cols.unplayed) = ~ismember(stimMeta.data{smCols.stimulus_id},...
%           prev_stim_ids);

  end
  
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
      if (isfield(params,'fun') && isfield(params.fun,'weight_fun') ...
              && isa(params.fun.weight_fun,'function_handle'))
          wfdata.sinfo = sinfo;
          wfdata.A = A;
          wfdata.W = W;
          wfdata.weight_cols = weight_cols;
          wfdata.stim_ids = stim_ids;
          wfdata.stimMeta = stimMeta;
          wfdata.smCols   = smCols;
          [A,W,Wi] = params.fun.weight_fun(wfdata,params);
      end
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

if (isfield(params,'fun') & isfield(params.fun,'pre_send') & ...
        isa(params.fun.pre_send,'function_handle'))
    psdata.A  = A;
    psdata.W  = W;
    psdata.prev_stim = prev_stim;
    try psdata.Wi = Wi;
    catch psdata.Wi = [];
    end
    [A,W,Wi,prev_stim.class]  = params.fun.pre_send(psdata,params);
elseif (isfield(Wi,'weights') & (length(Wi.weights) == 1))
    W = Wi.weights{1};
end

if (strcmp(W,'end'))
    stimloc = '';
    stim_id = '';
    stim_id_str = 'end';
else
    
    % Pass the attribute matrix A, which is an M x N matrix where,
    %    M = number of stimuli and 
    %    N = number of attributes
    % along with an attribute weighting vector, to select_from_weighted() which
    % draws randomly from the distribution of stimulus probabilities.
    item_idx = select_from_weighted(A,W);

    % Set the output variables
    stimloc = stimMeta.data{smCols.location}{item_idx};
    stim_id = stimMeta.data{smCols.stimulus_id}(item_idx);
    stim_id_str = num2str(stim_id);

    prev_stim.stimloc = stimloc;
    prev_stim.stim_id = stim_id;
    prev_stim.stim_id_str = stim_id_str;

    % Update any columns that need to be updated based on this selection
    if (isfield(params,'selection') && isfield(params.selection,'replacement')...
            && ~params.selection.replacement)
        A(item_idx,:) = 0;
    end
end

return
