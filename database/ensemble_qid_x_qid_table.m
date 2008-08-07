function an_st = ensemble_qid_x_qid_table(data_st,params)
% Produces a table comparing two question IDs.
% an_st = ensemble_qid_x_qid_table(data_st,params);
%
% Produces a table in which the levels for one question ID form the columns and
% the levels for the other question ID form the rows

% 05/18/07 PJ - added conn_id support, and resolution of {'display','report'} ambiguity

try
  conn_id = params.ensemble.conn_id;
catch
  conn_id = 1;
  tmp_conn_id = 1;
end

if isfield(params,'display')
  report_str = 'display';
else
  report_str = 'report';
end

an_st = {};
na = 0;

% Set the column constants
incol = set_var_col_const(data_st.vars);

% Make sure that we have either compqid or question and subquestion ID information
if ~isfield(incol,'compqid') & ~all(isfield(incol,{'question_id','subquestion'}))
  fprintf(['Did not find necessary question and subquestion or composite question ID' ...
	' information in the input data']);
  return
end

% Generate the compqid for each row in the data if necessary
if ~isfield(incol,'compqid')
  data_st.vars{end+1} = 'compqid';
  data_st.data{end+1} = make_compqid(data_st.data{incol.question_id}, ...
      data_st.data{incol.subquestion});
  incol = set_var_col_const(data_st.vars);
end

% Apply any specified filtering to the input data
if isfield(params,'filt')
  fprintf('Applying filtering criteria\n')
  data_st = ensemble_filter(data_st, params.filt);
end

% Make sure we only have two compqids
qids = unique(data_st.data{incol.compqid});
nqid = length(qids);

if nqid < 2
  fprintf('ensemble_qid_x_qid_table: Too few qids\n');
  return
elseif nqid > 2
  fprintf('ensemble_qid_x_qid_table: Too many (%d) qids\n', nqid);
  return
end

% Figure out how many categories we have in each of the questions
qinfo = mysql_extract_metadata('table','question',...
    'question_id',fix(qids), ...
    'conn_id', conn_id);
qinfo(~ismember([qinfo.compqid],qids)) = [];

rowqid_idx = find([qinfo.compqid] == params.table.qid_order(1));
colqid_idx = find([qinfo.compqid] == params.table.qid_order(2));
rowcats = qinfo(rowqid_idx).enum_values;
colcats = qinfo(colqid_idx).enum_values;
rowcats{end+1} = 'Total';
colcats{end+1} = 'Total';

nrows = length(rowcats);
ncols = length(colcats);

subids = unique(data_st.data{incol.subject_id});
nsubs = length(subids);

cntmtx = zeros(nrows,ncols,nsubs);
propmtx = zeros(nrows,ncols,nsubs);

% Initialize an analysis structure
na = na+1;

an_st{na} = init_analysis_struct;
an_st{na}.type = 'table_qid_x_qid_x_sub';
an_vars = {'subject_id','row_categories','col_categories','count','prop','nstim'};
an_st{na}.vars = an_vars;
an_cols = set_var_col_const(an_vars);
  
an_st{na}.data{an_cols.subject_id} = subids;  
an_st{na}.data{an_cols.row_categories} = rowcats;  
an_st{na}.data{an_cols.col_categories} = colcats;  

% Loop over subjects and stimuli to figure out how many items fall into each
% question combination
for isub = 1:nsubs
  filt.include.any.subject_id = subids(isub);
  subdata = ensemble_filter(data_st,filt);
  
  stimids = unique(subdata.data{incol.stimulus_id});
  nstim = length(stimids);
  
  % Generate the question ID masks
  rowqid_mask = subdata.data{incol.compqid} == params.table.qid_order(1);
  colqid_mask = subdata.data{incol.compqid} == params.table.qid_order(2);
  
  for istim = 1:nstim
    stim_mask = subdata.data{incol.stimulus_id} == stimids(istim);
    
    row_idx = ...
	enum2data(subdata.data{incol.response_enum}(stim_mask&rowqid_mask));
    col_idx = ...
	enum2data(subdata.data{incol.response_enum}(stim_mask&colqid_mask));
    
    cntmtx(row_idx,col_idx,isub) = cntmtx(row_idx,col_idx,isub)+1;
    
    % Update marginal counts
    cntmtx(row_idx,ncols,isub) = cntmtx(row_idx,ncols,isub)+1;
    cntmtx(nrows,col_idx,isub) = cntmtx(nrows,col_idx,isub)+1;
  end % for istim
  
  % Count the total number of stims
  cntmtx(nrows,ncols,isub) = sum(cntmtx(:,ncols,isub));
    propmtx(:,:,isub) = cntmtx(:,:,isub)/cntmtx(nrows,ncols,isub);
end % for isub=

an_st{na}.data{an_cols.count} = cntmtx;
an_st{na}.data{an_cols.prop} = propmtx;
an_st{na}.data{an_cols.nstim} = squeeze(cntmtx(nrows,ncols,:));
an_st{na}.report.printfun = @print_tables;

% Print the info if desired
if params.(report_str).tables.print | params.(report_str).tables.write2file
  an_st{na}.report.printfun(an_st{na},params.(report_str).tables);
end

if(exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end

end % ensemble_qid_x_qid_table

function an_st = init_analysis_struct
  an_st.type = '';
  an_st.vars = {};
  an_st.data = {};
  an_st.meta = [];
  
end % function an_st = init_analysis_struct

function print_tables(an_st,params)

  fid = 1;
  if params.write2file
    fid = fopen(params.fname,'wt');
    if fid ~= -1
      fprintf('Writing table to file: %s\n', params.fname);
    else
      fid = 1;
    end
  end

  an_cols = set_var_col_const(an_st.vars);
  colcats = an_st.data{an_cols.col_categories};
  rowcats = an_st.data{an_cols.row_categories};
  nrows = length(rowcats);
  
  nsub = length(an_st.data{an_cols.subject_id});
  
  fprintf(fid,'\n\nTotal number of instances\n');
  fprintf(fid,'\n%s\n', sprintf('\t%s', colcats{:}));
  for irow = 1:nrows
    data = sum(squeeze(an_st.data{an_cols.count}(irow,:,:)),2);
    fprintf(fid,'%s%s\n', rowcats{irow}, sprintf('\t%d',data));
  end
  
  fprintf(fid,'\n\nAverage percentage (mean +/- sem)\n');
  fprintf(fid,'\n%s\n', sprintf('\t%s', colcats{:}));
  for irow = 1:nrows
    data = mean(squeeze(an_st.data{an_cols.prop}(irow,:,:)),2);
    err = std(squeeze(an_st.data{an_cols.prop}(irow,:,:)),[],2)/sqrt(nsub);
    fprintf(fid,'%s%s\n', rowcats{irow}, sprintf('\t%2.1f (%2.1f)', [data err]'*100));
  end  
  
  if fid ~=1
    fclose(fid);
  end
end