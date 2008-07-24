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

% 09/14/05 Petr Janata


% Connect to host
try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
  tmp_conn_id = 1;
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

% Close the mysql connection if this was a temporary opening of the database
if exist('tmp_conn_id','var')
  mysql(conn_id,'close');
end
