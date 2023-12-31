function qtxt2qdf(host,database)
% Copies data from question.question_text to question_x_data_format.heading
% in the specified host and database
% qtxt2qdf(host,database);
%
% This function was created to contend with the issue that in a multipart
% question, the parent question text would not make it into the subheading
% field like the sub-questions.
%
% Copies data from question_text field in question table to heading field in
% the question_x_data_format table.
%
% NOTE: You must own this script to execute it. Use with caution

% 06/15/10 Petr Janata - turned script into function and added better
% connection handling

% Connect to host
params = mysql_login(struct('host',host,'database',database,'login_type','researcher'));
conn_id = mysql_make_conn(params);

% Get the question IDs and text 
mysql_str = sprintf('SELECT question_id, question_text FROM question;');
[qid,qt] = mysql(conn_id,mysql_str);

nq = length(qid);
for iq = 1:nq
  new_str = strrep(qt{iq},'"','\"');

  mysql_str = sprintf(['SELECT heading FROM question_x_data_format WHERE ' ...
	'question_id="%d" AND subquestion="1";'], qid(iq));
  curr_str = mysql(conn_id,mysql_str);
  if ~isempty(curr_str)
    curr_str = curr_str{1};
  else
    curr_str = '';
  end
  
  if isempty(curr_str)
    mysql_str = sprintf(['UPDATE question_x_data_format SET question_x_data_format.heading="%s" ' ...
	  'WHERE question_id="%d" AND subquestion="1";'], new_str, qid(iq));
    dummy = mysql(conn_id,mysql_str);
  else
    fprintf('Did not overwrite existing text for question (%d): %s\n', qid(iq), curr_str);
  end
end

% Close the mysql connection
mysql(conn_id,'close');
