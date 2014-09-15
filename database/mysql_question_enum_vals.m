function enum_vals = mysql_question_enum_vals(question_id, subquestion, conn_id)
% Returns a cell array of strings containing the enum labels for a given
% combination of question_id and subquestion

% 14Sep2014 Petr Janata
if nargin < 3
  error('%s: Requires 3 input arguments', mfilename)
end

if isempty(subquestion)
  subquestion = 1;
end

mysql_str = sprintf(['SELECT df.enum_values FROM data_format AS df, ' ...
  'question_x_data_format AS qdf ' ...
  'WHERE df.data_format_id = qdf.answer_format_id AND qdf.question_id = %d ' ...
  'AND qdf.subquestion = %d;'], question_id, subquestion);

enum_vals = mysql(conn_id, mysql_str);

enum_vals = regexp(regexprep(enum_vals{1},'"',''),',','split');

return