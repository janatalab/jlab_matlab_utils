function an_st = ensemble_data_by_question(data_st,params)
% Filters/reshapes a struct such that compqids are now fields within the return struct.
% an_st = ensemble_data_by_question(data_st,params);
%
% Filters and reshapes a data structure such that desired composite question
% IDs are now the different variables in the resulting data structure.
%
% params.compqids - list of composite qids that we want to extract. If none are
%                   given, then use all.
% params.extract_vars - optional cell string array which gives list of
%                       variables in data_st that should be transferred to the
%                       new structure. Default is to transfer all except for
%                       compqid

% 02/04/07 Petr Janata
% 10/07/08 PJ - generalized to handle other databases
% 22/07/10 PJ - fixed conn_id handling

an_st = ensemble_init_data_struct;
an_st.type = 'data_by_compqid'; 

% Set the column constants
incol = set_var_col_const(data_st.vars);

% Make sure that we have either compqid or question and subquestion ID information
if ~isfield(incol,'compqid') && ~all(isfield(incol,{'question_id','subquestion'}))
  fprintf(['Did not find necessary question and subquestion or composite question ID' ...
	' information in the input data']);
  return
end

% Extract info regarding the database and connection ID we should be talking to
param_fld_names = {'ensemble','mysql'};
idxs = find(isfield(params,param_fld_names));
if isempty(idxs)
  error('%s: Do not have sufficient database connection information', mfilename)
else
  database = params.(param_fld_names{idxs(1)}).database;
  conn_id = params.(param_fld_names{idxs(1)}).conn_id;
end

% Apply any specified filtering to the input data
if isfield(params,'filt')
  fprintf('Applying filtering criteria\n')
  data_st = ensemble_filter(data_st, params.filt);
end

% Generate the compqid for each row in the data if necessary.  We want to do
% this before filtering, since some of the filters might be specified in terms
% of the composite question ID (compqid).  Because we need these anyway, we
% might as well make sure we have them now.
if ~isfield(incol,'compqid')
  data_st.vars{end+1} = 'compqid';
  data_st.data{end+1} = make_compqid(data_st.data{incol.question_id}, ...
      data_st.data{incol.subquestion});
  incol = set_var_col_const(data_st.vars);
end

% Check to see if we need to prune the dataset further to contain only the
% composite question IDs that we want.
unique_src_compqid = unique(data_st.data{incol.compqid}); % existing compqids
if isfield(params,'compqids') && ~isempty(params.compqids)
  if ~(all(ismember(unique_src_compqid,params.compqids)))
    filt = [];  % clear the filter structure
    filt.include.all.compqid = params.compqids;
    data_st = ensemble_filter(data_st,filt);
  end
  compqids = params.compqids;
else
  compqids = unique(data_st.data{incol.compqid});
end

% Get final list of compqids
nqid = length(compqids);

% Extract the question info and attach it to the output metadata
qinfo = mysql_extract_metadata('table','question', ...
    'question_id',unique(fix(compqids)),'conn_id',conn_id);
qinfo = qinfo(ismember([qinfo.compqid],compqids));
an_st.meta.question = qinfo;

% Make list of variables to extract 
if isfield(params,'extract_vars') && ~isempty(params.extract_vars)
  extract_vars = params.extract_vars;
else
  extract_vars = data_st.vars;
end

% get column indices into the source data
[dummy,src_col] = ismember(extract_vars, data_st.vars);

for iqid = 1:nqid
  an_st.vars{iqid} = sprintf('compqid_%d_%d',qinfo(iqid).question_id,qinfo(iqid).subquestion);

  % Initialize the data structure that will be assigned to this variable
  curr_st = ensemble_init_data_struct;
  curr_st.type = 'raw_data';
  curr_st.vars = extract_vars;
  curr_st_cols = set_var_col_const(curr_st.vars);
    
  % Filter the data for this particular compqid
  clear filt
  filt.include.all.compqid = qinfo(iqid).compqid;
  tmpdata = ensemble_filter(data_st,filt);

  % Retain only the source columns we want
  curr_st.data = tmpdata.data(src_col);  

  an_st.data{iqid} = curr_st;
end % for iqid

