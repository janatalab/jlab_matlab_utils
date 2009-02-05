function [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim_v8(varargin)

% [stim_id_str,stim_id,stimloc] = get_ensemble_weighted_stim_v8(varargin);
%
% Returns the ID of a stimulus to present as a string.
%
% This is intended to be a very generic script that to manage stimulus
% selection needs of Ensemble experiments.  Experiment specific parameters
% should be stored in a paramters file, e.g. vdc_prescreen_params_v1.m, the
% name of which is passed into the function following a params_file tag.
%
% INPUT
%  'reset' - this will clear all persistent variables and return
%    immediately. it is useful for debugging or simulation circumstances
%    when one would want to run multiple subjects within one matlab session
% 
%  subid - subject_id field for response_table
%
%  params.mysql - must contain host,database,conn_id fields with str,str,int values
%  params.resptbl_name - (str) name of this experiment's ensemble database response table
%  params.version - (str) param file version info to be input to misc_info col in resptbl_name
%  params.targ_form_id - (int) form_id whose row will have version input into the misc_info
%  params.prev_expts - (cell array of strings) stim presented during these
%    experiments will be excluded from the current experiment.
%  params.return_stim_during_init - (int) if set to 1, this script will begin
%    returning stim upon first call
%
%  params.fun.stim_compiler - (fh) default is mysql_get_stim_by_attribute, see
%    this for an example of expected final output.
%  params.fun.weight_fun - (fh) only executed if params.selection.type = 'weighted'
%    this is a function handle to an experiment-specific function that handles
%    apriori stim weighting, through generation of A (stim*attribute weight
%    matrix), W (attribute weight vector to be applied to A) and Wi (a struct
%    containing an array of possible weight vectors in Wi.weights and an order
%    for those weights to be applied in Wi.order ... contents of Wi.order
%    should be the index of the Wi.weights vector to be used to select stim at
%    any given index of Wi.order. Indices of Wi.order are equivalent to the
%    number of stim already presented + 1. In the case that Wi is generated and
%    used, the first value of W shall be Wi.weights{Wi.order{1}}. Wi.counter
%    should maintain the current index of Wi.order that is being used.). If you
%    are using Wi, W and Wi should interact through params.fun.pre_send.
%  params.fun.pre_send - (fh) pre_send is a function handle to an
%    experiment-specific funtion to carry out any calculations or database
%    operations necessary before choosing and returning the next stimulus to be
%    presented. This should be used particularly when you have created Wi and
%    wish to select the next Wi.weights and return it in the form of W.
%
%  params.weights - (struct) fieldnames should all be names of weights, and values of
%    the initial value of these weights for all stimuli.
%  params.weights.init_prob - the intial probability assigned to all stimuli,
%    if params.selection.type = 'equal_probability'
%
%  params.selection.replacement - if set to 0 (no replacement), all previous
%    stim found by params.prev_expts will be given a 0 probability to be
%    played, and each stim presented will be given a 0 probability for future
%    stimulus selection.
%  params.selection.type - (str) legal values are 'weighted' and 'equal_probability'
%

% Modification history
% 
% FB Sep 27 2007 - adding functionality to weight stims, as well as to 
% dynamically update stim weighting ... this has been moved out of the
% 'is_initialized' conditional, so that it can be continuously updated
% 
% FB Feb 04 2008 - documentation, cleaning up code
%
% FB Apr 21 2008 - added 'reset' function
%
% 12/27/08 v7 PJ - exclusion of previously encountered stimuli was not being
%                  properly handled.  It should take place after the stimulus
%                  weights have been set.
% 
% 2009.02.04 FB v8 - when dealing with function handles for stim_compiler,
% weight_fun or pre_send, if the provided var is a string, use str2func,
% but if it is a function_handle, then treat it as such

script_version = 'get_ensemble_weighted_stim_v6.0';

stim_id_str = '';
stim_id = '';
stimloc = '';

% Parse input parameters
for iarg = 1:2:nargin
  curr_arg = varargin{iarg};
  
  switch curr_arg
      case {'params','params_file','param_file'}
          param_file = varargin{iarg+1};
      case {'subid','subject_id','subjectid'}
          subid = varargin{iarg+1};
      case {'reset'}
          clear is_initialized conn_id params stimids stimMeta smCols
          clear A W Wi weight_cols prev_stim
          return
    
  end % swich curr_arg
end % for iarg

if ~exist('param_file') || isempty(param_file)
  fprintf('No parameter file specified\n');
  return
end

if ~exist('subid')
  fprintf('No subid specified\n')
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
  
  try debug = params.debug; catch debug = false; end

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
  %
  % NOTE: This is a dangerous way of doing things because it is possible for
  % multiple experiments to write to the same response table. It's not clear by
  % the experiment isn't specified by name. PJ.
  mysql_str = sprintf('SELECT experiment_id FROM experiment WHERE response_table="%s";',...
      params.resptbl_name);
  curr_expid = mysql(conn_id, mysql_str);

  %   mysql(conn_id,sprintf('use %s',params.mysql.database));
  % Insert the script name into the table
  if (isstruct(debug) && isfield(debug,'on') && debug.on ~= false) || (isstruct(debug) ...
	&& ~isfield(debug,'on')) || (~isstruct(debug) && ~debug)
    insert_str = sprintf(['UPDATE %s SET %s.misc_info = ''param_file:%s'' ' ...
	  'WHERE subject_id=''%s'' AND experiment_id=%d AND form_id=%d;'], ...
	params.resptbl_name,params.resptbl_name,params.version,...
	subid,curr_expid,params.targ_form_id);
  
    mysql(conn_id,insert_str);
  end

  % Get info on this subject including the experiments they've been in
  subinfo = mysql_get_sinfo(subid,conn_id);
  subinfo.subid = subid;

  % if stim_compiler is defined, then run that function, else get the
  % attributes as they are defined within attrib_names, through
  % mysql_get_stim_by_attribute
  if (isfield(params,'fun') & isfield(params.fun,'stim_compiler'))
      if ischar(params.fun.stim_compiler)
        scfh = str2func(params.fun.stim_compiler);
      elseif iscell(params.fun.stim_compiler)
        if length(params.fun.stim_compiler) ~= 2 || ...
                ~ischar(params.fun.stim_compiler{1}) || ...
                ~ischar(params.fun.stim_compiler{2})
          error(['must provide library name and sub-function for '...
              'stim compiler function extraction in string form\n']);
        end
        sclfh = str2func(params.fun.stim_compiler{1});
        scfh = sclfh(params.fun.stim_compiler{2});
      elseif isa(params.fun.stim_compiler,'function_handle')
        scfh = params.fun.stim_compiler;
      else
        error(['params.fun.stim_compiler must either be a function '...
            'handle or a string-name of a valid function\n']);
      end

      [stimdata,stimvars] = scfh(params);
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
      if (isfield(params,'fun') && isfield(params.fun,'weight_fun'))
          if ischar(params.fun.weight_fun)
            wffh = str2func(params.fun.weight_fun);
          elseif iscell(params.fun.weight_fun)
            if length(params.fun.weight_fun) ~= 2 || ...
                    ~ischar(params.fun.weight_fun{1}) || ...
                    ~ischar(params.fun.weight_fun{2})
              error(['must provide library name and sub-function for '...
                  'weight fun function extraction in string form\n']);
            end
            wflfh = str2func(params.fun.weight_fun{1});
            wffh = wflfh(params.fun.weight_fun{2});
          elseif isa(params.fun.weight_fun,'function_handle')
            wffh = params.fun.weight_fun;
          else
            error(['params.fun.weight_fun must either be a function '...
                'handle or a string-name of a valid function\n']);
          end

          wfdata.sinfo = subinfo;
          wfdata.A = A;
          wfdata.W = W;
          wfdata.weight_cols = weight_cols;
          wfdata.stim_ids = stim_ids;
          wfdata.stimMeta = stimMeta;
          wfdata.smCols   = smCols;
          [A,W,Wi] = wffh(wfdata,params);
      end
  end

  if (isfield(params,'selection') && isfield(params.selection,'replacement') && ...
          ~params.selection.replacement)

      % % % %   
      % get stim this participant has already heard

      % Figure out if the subject was in previous versions of this
      % experiment, and get a list of stimuli that were presented to the
      % subject
      try prev_expts = params.prev_expts; catch prev_expts = ''; end
      if ~isempty(prev_expts)
	prev_exp_idxs = find(ismember(subinfo.exp_names, params.prev_expts));
	prev_stim_ids = find_presented_stims(subinfo, prev_exp_idxs, conn_id);
	prev_stim_idxs = find(ismember(stimMeta.data{smCols.stimulus_id},...
	    prev_stim_ids));
	A(prev_stim_idxs,:) = 0;
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

if (isfield(params,'fun') & isfield(params.fun,'pre_send'))
  if ischar(params.fun.pre_send)
    psfh = str2func(params.fun.pre_send);
  elseif iscell(params.fun.pre_send)
    if length(params.fun.pre_send) ~= 2 || ...
            ~ischar(params.fun.pre_send{1}) || ...
            ~ischar(params.fun.pre_send{2})
      error(['must provide library name and sub-function for '...
          'pre-send function extraction in string form\n']);
    end
    pslfh = str2func(params.fun.pre_send{1});
    psfh = pslfh(params.fun.pre_send{2});
  elseif isa(params.fun.pre_send,'function_handle')
    psfh = params.fun.pre_send;
  else
    error(['params.fun.pre_send must either be a function '...
        'handle or a string-name of a valid function\n']);
  end

  psdata.A  = A;
  psdata.W  = W;
  psdata.prev_stim = prev_stim;
  try psdata.Wi = Wi;
  catch psdata.Wi = [];
  end
  [A,W,Wi,prev_stim.class]  = psfh(psdata,params);
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

    if isempty(item_idx)
        stimloc = '';
        stim_id = '';
        stim_id_str = 'end';
    else
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
end

return
