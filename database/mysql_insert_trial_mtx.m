function result = mysql_insert_trial_mtx(trial_mtx,params)
% Inserts trial information into the trial table.
%
% USAGE: result = mysql_insert_trial_mtx(trial_mtx,params);
%
% trial_mtx: stimulus_id1, stimulus_id2, question_id, data_format_id, correct_response_enum

% 07/18/09 Petr Janata

% Check for connection to database
try 
  conn_id = params.conn_id;
  tmp_conn_id = 0;
catch
  tmp_conn_id = 1;
  conn_id = 0;
end

if mysql_check_conn(conn_id)
  try host = params.host; catch host = []; end
  try database = params.database; catch database = []; end
  mysql_make_conn(host,database,conn_id);
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

mysql_str = sprintf(['SELECT trial_id, stimulus_id1, stimulus_id2 question_id data_format_id '  ...
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
elseif length(trial_ids) < num_trial_types
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
  mysql(conn_id,'close');
end
