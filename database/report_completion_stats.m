function rp = report_completion_stats(einfo, rp)
% Generates a report of completion statistics for given experiment information.
% rp = report_completion_stats(expinfo, report_params)
%
% Generates a report on the completion statistics for subjects associated with
% an experiment described in the expinfo structure.  
%
% For each subject (and each session), the session information is consulted in
% the response table to determine the last form submitted in that session.
%
% The rp structure can contain several fields that guide the behavior of this
% analysis module.
%
% .ticket_code -- if this cell array of strings is present, only those subjects
%                 who accessed the experiment with a ticket from this list are
%                 included in the report
%
% .print_report -- should the report be printed to the screen
%
% See also: mysql_get_expinfo()

% 10/29/06 Petr Janata - started script

% Check for an active mysql connection and open one if necessary
try rp.conn_id(1);
  conn_id = rp.conn_id;
catch
  conn_id = 0;
end
mysql_check_conn(conn_id,'open');

nexp = length(einfo);
if nexp > 1
  error(sprintf('Too many (%d) experiments\n', nexp))
end

% Check to see if we are restricting the search by ticket
if isfield(rp,'ticket_code') && ~isempty(rp.ticket_code)
  ticket_code = rp.ticket_code;
  if ~iscell(ticket_code)
    ticket_code = {ticket_code};
  end
  ntickets = length(ticket_code);

  ticket_str = sprintf('"%s"', cell2str(ticket_code,'","'));
  mysql_str = sprintf(['SELECT session.session_id FROM session, ticket WHERE' ...
	' session.ticket_id = ticket.ticket_id AND ticket.ticket_code IN (%s);'], ticket_str);
  [ticket_sessids] = mysql(conn_id,mysql_str);
else
  ticket_sessids = [];
end

if ~isfield(rp,'report_str')
  rp.report_str = {};
end

fid = 1;
if isfield(rp,'report_fname')
  fid = fopen(rp.report_fname,'wt');
  if fid == -1
    fprintf('Problem opening reporting file: %s\n', rp.report_fname);
    fid = 1;
  end
end
    

if isfield(rp,'print_report') && rp.print_report
  PRINT_REPORT = 1;
  % Print the header
  fprintf(fid,['Subject Name\tSubject ID\tSession ID\tStart Time\tElapsed time' ...
	' (h)\tLast Form ID\tLast Form Name\n']);
else
  PRINT_REPORT = 0;
end

cexp = einfo(1);  % make a copy of the current experiment info

% Get the experiment form list
exp_formid_list = cexp.forms.ids;
exp_formname_list = cexp.forms.names;

% Set up a vector to tally the number of times a form was a terminal form
terminal_form_count = zeros(size(exp_formid_list));

nsubs = length(cexp.subs.names);
for isub = 1:nsubs
  % Determine the number of sessions/subject
  nsess = cexp.subs.nsess(isub);
  
  for isess = 1:nsess
    sess_id = cexp.subs.sess{isub}(isess);
    
    if ~isempty(ticket_sessids) && ~any(ticket_sessids == sess_id)
      continue
    end
    
    % Get the session time stamp
    mysql_str = sprintf(['SELECT date_time,end_datetime FROM session ' ...
	  'WHERE session_id=%d;'], sess_id);
    [start_datenum, end_datenum] = mysql(conn_id, mysql_str);
    
    % Get a list of forms that the subject completed
    mysql_str = sprintf(['SELECT DISTINCT form_id FROM %s '...
	  'WHERE session_id=%d;'], cexp.resp_table, sess_id);
    completed_ids = mysql(conn_id,mysql_str);
    
    % See how far the subject got
    last_form_idx = max(find(ismember(exp_formid_list,completed_ids)));
    
    rp.subs.last_form_idx{isub}(isess) = last_form_idx;
    rp.subs.last_form_id{isub}(isess) = exp_formid_list(last_form_idx);
    rp.subs.last_form_name{isub}{isess} = exp_formname_list{last_form_idx};
    rp.subs.starttime{isub}(isess) = start_datenum;
    if ~isempty(end_datenum)
      elapsed_time = etime(datevec(end_datenum),datevec(start_datenum));
      rp.subs.elapsedtime{isub}(isess) = elapsed_time;
    else
      elapsed_time = -1;
    end
    
    % Increment the tally
    terminal_form_count(last_form_idx) = terminal_form_count(last_form_idx)+1;

    curr_str = sprintf('%20s\t%s\t%d\t%s\t%1.2f\t%d\t%s\n', ...
	cexp.subs.names{isub}, ...
	cexp.subs.ids{isub}, ...
	sess_id, ...
	datestr(start_datenum), ...
	elapsed_time/3600, ...
	exp_formid_list(last_form_idx), ...
	exp_formname_list{last_form_idx});
    rp.report_str{end+1} = curr_str;
    
    if PRINT_REPORT
      fprintf(fid,'%s', curr_str);
    end
  end % for isess
end % for isub=

rp.terminal_form_count = terminal_form_count;

if ~conn_id
  mysql(conn_id,'close')
end

% Print the stuff out if desired
