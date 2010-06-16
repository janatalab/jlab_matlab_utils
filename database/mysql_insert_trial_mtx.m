function result = mysql_insert_trial_mtx(trial_mtx,params)
% Inserts trial information into the trial table.
%
% USAGE: result = mysql_insert_trial_mtx(trial_mtx,params);
%
% trial_mtx: stimulus_id1, stimulus_id2, question_id, data_format_id, correct_response_enum
%
% params should be a structure containing database information so that it can connect to the
% database: 
% Either: host,database,user,passwd
% or conn_id

% 07/18/09 Petr Janata
% 06/15/10 PJ - sanitized mysql_make_conn

result = [];

if mysql_check_conn(params.conn_id)
  if ~all(isfield(params,{'host','database'}))
    error('%s: Insufficient information to establish database connection', mfilename);
  else
    params.login_type = 'researcher';
    params = mysql_login(params);
  end
  params = mysql_make_conn(params);
  tmp_conn_id = 1;
else
  tmp_conn_id = 0;
end

% 
% Check to see if trial IDs have already been created in the trial table
%
stim1_str = sprintf('%d,', unique(trial_mtx(:,1)));
stim1_str(end) = '';

stim2_str = sprintf('%d,', unique(trial_mtx(:,2)));
stim2_str(end) = '';

qid_str = sprintf('%d,', unique(trial_mtx(:,3)));
qid_str(end) = '';

dfid_str = sprintf('%d,', unique(trial_mtx(:,4)));
dfid_str(end) = '';

mysql_str = sprintf(['SELECT trial_id, stimulus_id1, stimulus_id2, question_id, data_format_id '  ...
      'FROM trial WHERE ' ...
      'stimulus_id1 IN (%s) AND ' ...
      'stimulus_id2 IN (%s) AND ' ...
      'question_id IN (%s) AND ' ...
      'data_format_id IN (%s);'],...
    stim1_str, stim2_str, qid_str, dfid_str);
[trial_ids, stim1_ids, stim2_ids, question_ids, dfids] = mysql(params.conn_id,mysql_str);

%
% Need to check if any trial IDs need to be inserted
%
insert_mtx = [];
if isempty(trial_ids)
  insert_mtx = trial_mtx;
elseif length(trial_ids) < size(trial_mtx,1)
  % Remove existing trials
  for itr = 1:length(trial_ids)
    existing_idxs = ...
	trial_mtx(:,1)==stim1_ids(itr) & ...
	trial_mtx(:,2)==stim2_ids(itr) & ...
	trial_mtx(:,3)==question_ids(itr) & ...
	trial_mtx(:,4)==dfids(itr);
    
    trial_mtx(existing_idxs,:) = [];
  end
  insert_mtx = trial_mtx;
end

if ~isempty(insert_mtx)
  fprintf('Inserting %d trials into trial table\n', size(insert_mtx,1));
  value_str = sprintf('("%d","%d","%d","%d","%d"),', insert_mtx');
  value_str(end) = '';

  mysql_str = sprintf(['INSERT INTO trial (' ...
	'stimulus_id1, ' ...
	'stimulus_id2, ' ...
	'question_id, ' ...
	'data_format_id, ' ...
	'correct_response_enum) ' ...
	'VALUES %s;'], value_str);
  
  result = mysql(params.conn_id, mysql_str);
end

% Close connection if necessary
if tmp_conn_id
  mysql(params.conn_id,'close');
end
