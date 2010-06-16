function [dfid] = mysql_qid2dfid(qid, conn_id)
% Returns dfid associated with given qid
%
% [dfid] = mysql_qid2dfid(qid);
%
% Returns the data format id (dfid) associated with the unique
% question/subquestion combinations provided in qid
% 
% Column 1 of qid contains questions, and column 2 contains the associated
% subquestions.  A single part question has a subquestion id of 1
%
% conn_id - connection to database - required

% 09/14/05 Petr Janata
% 06/15/10 PJ - mysql_make_conn sanitization

% Check for valid connection to database
if ~exist('conn_id','var') || isempty(conn_id) || mysql(conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

[unique_quest_ids, quest_idxs] = unique(qid,'rows'); % get the unique questions
nquest = size(unique_quest_ids,1);

% Create the mysql query string
qid_str = sprintf('(question_id=%d AND subquestion=%d) OR ', unique_quest_ids');
qid_str(end-3:end) = [];  

mysql_str = sprintf(['SELECT type, data_format_id FROM data_format ' ...
      'RIGHT JOIN question_x_data_format ON' ...
      ' data_format.data_format_id=question_x_data_format.answer_format_id ' ...
      'WHERE (%s);'], qid_str);
[types, dfid] = mysql(conn_id,mysql_str);

return