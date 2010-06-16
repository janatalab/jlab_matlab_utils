function [subids, sessids, form_id, question_id, timestamps] = check_completion(tbl_name,varargin)
% Returns a list of subject and session IDs that responded to a given form or question.
%
% [subids, sessids, form_id, question_id, timestamp] = check_completion(tbl_name,varargin);
%
% If no output arguments are specified, the results are
% printed to standard output.
%
% The form ID and/or question ID are specified by tag/value pairs, e.g.
% 
% 'form_id', 123
% 'question_id', 456

% 07/02/06 Petr Janata
% 06/15/10 PJ - eliminated mysql_make_conn

% Parse input arguments
narg = length(varargin);

cond_str = '';
FIND_BY_SUBJECT = 0;
FIND_BY_TICKET = 0;

conn_id = [];
for iarg = 1:2:narg
  curr_arg = varargin{iarg};
  if isempty(cond_str)
    prefix_str = '';
  else
    prefix_str = 'AND';
  end
  
  if isstruct(curr_arg)
    ap = curr_arg;
    if isfield(ap,'ticket_code') && ~isempty(ap.ticket_code)
      ticket_code = ap.ticket_code;
      if ~iscell(ticket_code)
        ticket_code = {ticket_code};
      end
      ntickets = length(ticket_code);
      FIND_BY_TICKET = 1;
    end
  else
    switch curr_arg
      case 'form_id'
        form_id = varargin{iarg+1};
        cond_str = sprintf('%s %s form_id=%d ', cond_str, prefix_str, form_id);
      case 'question_id'
        question_id = varargin{iarg+1}
        cond_str = sprintf('%s %s question_id=%d ', cond_str, prefix_str, question_id);
      case 'subject_id'
        subids = varargin{iarg+1};
        cond_str = sprintf('%s %s subject_id="%s" ', cond_str, prefix_str, subids);
        FIND_BY_SUBJECT = 1;
      case 'conn_id'
        conn_id = varargin{iarg+1};
    end
  end % if isstruct(curr_arg)
end % for iarg=

% Connect to host with a temporary connection if necessary
if isempty(conn_id) || mysql(conn_id,'status')
  error('%s: No connection ID specifed or connection not open')
end

% 
% If we are searching by ticket codes, get a listing of the appropriate sessions
% based on the list of ticket codes
%
if FIND_BY_TICKET
  ticket_str = sprintf('"%s"', cell2str(ticket_code,','));
  mysql_str = sprintf(['SELECT session.session_id FROM session, ticket WHERE' ...
	' session.ticket_id = ticket.ticket_id AND ticket.ticket_code IN (%s);'], ticket_str);
  [sessids] = mysql(conn_id,mysql_str);
end

%
% If we are searching by subject, return the form and question ID last completed
%
if FIND_BY_SUBJECT
  var_list = 'session_id, form_id, question_id, date_time';
else
  var_list = 'session_id, subject_id, date_time';
end
mysql_str = sprintf('SELECT %s FROM %s WHERE %s;', var_list, tbl_name, cond_str);

if FIND_BY_SUBJECT
  [sessids, formids, qids, timestamps] = mysql(conn_id, mysql_str);
  form_id = formids(end);
  question_id = qids(end);
  timestamps = timestamps(end);
  sessids = sessids(end);
else
  [sessids, subids, timestamps] = mysql(conn_id, mysql_str);
end

%
% Close the mysql connection if this was a temporary opening of the database
%

if ~conn_id
  mysql(conn_id,'close');
end

if ~nargout
  fprintf('SubjectID\tSession\tFormID\tQuestionID\tTimestamp\n');
  for is = 1:length(subids)
    sub = subids{is};
    sess = sessids(is);
    tstmp = timestamps(is);
    fprintf('%s\t%d\t%d\t%d\t%s\n', sub, sess, form_id, question_id, datestr(tstmp));
  end
end
