function ensemble_export_stims(data_st,params)
% Copies stimuli that are registered in the Ensemble database to a folder
%
% Stimuli to be exported can be specified either in stimulus_id variable in
% data_st or can be harvested from a response table associated with the
% experiment specified in params.experiment_name
%

% 23Feb2013 Petr Janata

% Establish a connection to the database
params.mysql.conn_id = mysql_make_conn(mysql_login(params));
conn_id = params.mysql.conn_id;

% Determine whether stimuli are specified in data_st or whether we need to
% retrieve them via the response table
if ~isempty(data_st) && ismember(data_st.vars,'stimulus_id')
  cols = set_var_col_const(data_st.vars);
  stimids = unique(data_st.data{cols.stimulus_id});
else
  srchFlds = {'expname','experiment_name','experiment_title'};
  expNameMask = isfield(params.ensemble,srchFlds);
  if ~any(expNameMask)
    error('Cannot locate experiment name in params structure');
  else
    experiment_title = params.ensemble.(srchFlds{find(expNameMask,1,'first')});
  end

  % Get the experiment_id
  mysql_str = sprintf(['SELECT experiment_id, response_table FROM experiment WHERE ' ...
    'experiment_title = "%s";'], experiment_title);
  [expID, respTable] = mysql(conn_id, mysql_str);
  
  % Retrieve unique stims associated with this response table and this
  % experiment I
  mysql_str = sprintf(['SELECT DISTINCT stimulus_id FROM %s ' ...
    'WHERE experiment_id = %d AND stimulus_id IS NOT NULL;'], respTable{1}, expID);
  stimids = mysql(conn_id, mysql_str);
  
end

% Now retrieve stimulus location info
fprintf('Retrieving information for %d stimuli\n', length(stimids));
stim_st = ensemble_get_stiminfo(struct('vars','stimulus_id','data',{{stimids}}),params);
scols = set_var_col_const(stim_st.vars);

% Now call the stim processing wrapper
params.funcName = 'ensemble_copy_stimulus';
params.funcParams.stimroot = params.ensemble.stimroot;
params.funcParams.outpath = params.paths.stimpath;

stim_st = stim_st.data{scols.stimulus_metadata};
ensemble_processStims(stim_st,params);

return
end