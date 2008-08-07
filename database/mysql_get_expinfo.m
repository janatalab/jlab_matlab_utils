function [expinfo] = mysql_get_expinfo(expmt, host, database, conn_id)
% Gets information for an experiment from a mysql database
%
% [expinfo] = mysql_get_expinfo(expmt, host, db, conn_id);
%
% expmt -- experiment to look for.  If this is left empty, a list of
% experiments is displayed to the screen
% host -- machine running the mysql server
% db -- database on the server
%

% 01/25/05 Petr Janata
% 06/26/05 PJ - modified conn_id handling
% 03/13/06 PJ - Exclusion of tmp_ entries when retrieving subject_id
% 06/20/06 PJ - Added retrieval of form names and form IDs
% 07/17/06 PJ - Accelerated script using GROUP BY clauses
% 10/28/06 PJ - Made connecting to arbitrary desired database a bit more
%               robust. Now returns information for list of experiments in
%               expmt.

expinfo = [];

% Do some parameter checking
msg = nargchk(0,4,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

if nargin < 4
  conn_id = 0;
end

% See if we have a connection
if mysql(conn_id,'status')
  try 
    mysql_make_conn(host, database, conn_id);
  catch  
    tmp_conn_id = 1;
    mysql_make_conn;
    conn_id = 0;
  end
end

% Get a list of experiment titles, IDs, and response tables
mysql_str = sprintf('SELECT experiment_id,experiment_title,response_table FROM experiment');
[exp_ids, exp_names, resp_table_list] = mysql(conn_id,mysql_str);

% Display experiment list if no experiment name passed in. Otherwise return
% info on desired experiment
if ~exist('expmt') || isempty(expmt)
  fprintf('%s\n',cell2str(exp_names,'\n'));
else
  if ~iscell(expmt)
    expmt = {expmt};
  end
  
  nexp = length(expmt);
  
  for iexp = 1:nexp
    expinfo(iexp).name = expmt{iexp};
    
    % Determine which row of the response corresponds to the desired experiment
    row_idx = find(ismember(exp_names,expmt{iexp}));
    
    expinfo(iexp).expid = exp_ids(row_idx);
    
    % Get the associated response table name
    resp_table = resp_table_list{row_idx};
    expinfo(iexp).resp_table = resp_table;

    % Get a list of the form IDs in this experiment
    sql_str = sprintf(['SELECT form_id, form_order FROM experiment_x_form WHERE' ...
	  ' experiment_id = %d GROUP BY form_id;'], expinfo(iexp).expid);
    [form_ids, form_order] = mysql(conn_id,sql_str);
    
    % Sort form_ids by form_orders
    [sorted_form_order, form_order_idxs] = sort(form_order);
    expinfo(iexp).forms.ids = form_ids(form_order_idxs);
    
    form_id_str = sprintf('%d,', form_ids(form_order_idxs));
    form_id_str(end) = [];
    
    sql_str = sprintf('SELECT form_name, form_id FROM form WHERE form_id IN (%s);', form_id_str);
    [form_names, form_ids2] = mysql(conn_id,sql_str);

    % Need to sort form_ids such that they line up, given that they were
    % retrieved from separate tables. This could probably be accomplished with
    % a single nice SQL query, but I don't have the manual in front of me right
    % now, so I'll do it a complicated Matlab way.  PJ
    
    [in,ia,ib] = intersect(form_ids(form_order_idxs),form_ids2);
    [sort_ia, sort_ia_idx] = sort(ia); % get original list lined back up
    form_order_idxs2 = ib(sort_ia_idx); % use to dereference 2nd list index
                                        % generated by intersect
    
    expinfo(iexp).forms.names = form_names(form_order_idxs2);
    
    % Get desired variables from the response table, but ignore entries beginning
    % with tmp_ which denote temporary subject IDs left over from aborted
    % consenting processes.
    sql_str = sprintf(['SELECT subject_id, session_id, ' ...
	  'MIN(date_time) AS start_time, MAX(date_time) AS stop_time ' ...
	  'FROM %s WHERE subject_id NOT LIKE "tmp_%%" GROUP BY session_id;'], resp_table);
    
    [ids,sess_ids,start_times,stop_times] = mysql(conn_id,sql_str);
%     start_times = datenum(start_times,'yyyy-mm-dd HH:MM:SS');
%     stop_times = datenum(stop_times,'yyyy-mm-dd HH:MM:SS');
    
            % Get a list of subjects in the experiment
    [expinfo(iexp).subs.ids] = unique(ids);
    nsubs = length(expinfo(iexp).subs.ids);
    
    % Get subject info
    sub_str = sprintf('''%s'',', expinfo(iexp).subs.ids{:});
    sub_str(end) = [];
    mysql_str = sprintf(['SELECT name_first, name_last, dob FROM subject WHERE' ...
	  ' subject_id IN (%s)'], sub_str);
    [first_names, last_names, birthdays] = mysql(conn_id, mysql_str);
    
    % Determine the number of sessions for each subject in the experiment
    for isub = 1:nsubs
      sub_idx_mask = ismember(ids, expinfo(iexp).subs.ids{isub});
      [curr_sess_ids, curr_idxs] = unique(sess_ids(sub_idx_mask));
      expinfo(iexp).subs.sess{isub} = curr_sess_ids;
      expinfo(iexp).subs.nsess(isub) = length(curr_sess_ids);
      expinfo(iexp).subs.names{isub} = sprintf('%s %s',first_names{isub},last_names{isub});
      datediff = datevec(start_times(curr_idxs)-birthdays(isub));

      % take the average age if dealing with multiple sessions 
      expinfo(iexp).subs.ages(isub) = mean(datediff(:,1));  
    end % for isub
    
    % Deal with things from the perspective of sessions
    expinfo(iexp).sess.ids = sess_ids;
    expinfo(iexp).sess.subids = ids;    
    expinfo(iexp).sess.start_time = start_times;
    expinfo(iexp).sess.stop_time = stop_times;
    
    % Deal with the rare situation (from early days of Ensemble) where we don't
    % have a session stop time.
    missing_idxs = find(isnan(expinfo(iexp).sess.stop_time));
    if any(missing_idxs)
      sess_str = sprintf('%d,',expinfo(iexp).sess.ids(missing_idxs));
      sess_str(end) = [];
      mysql_str = sprintf(['SELECT session_id, MAX(date_time) AS endtime FROM %s ' ...
	    'WHERE session_id IN (%s) GROUP BY session_id'], resp_table, sess_str);
      tic
      [sess_ids, end_times] = mysql(conn_id, mysql_str);
      toc
      expinfo(iexp).sess.stop_time(missing_idxs) = end_times;
    end % if any(missing_idxs)
  end % for iexp
end % if ~exist('expmt') || isempty(expmt)

% Close the mysql connection if this was a temporary opening of the database
if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end