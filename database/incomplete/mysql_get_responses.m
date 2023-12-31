function responses = mysql_get_responses(subject_id,response_table)
%
% Join the form, question, and data format records to all of the responses from a given subject in a response table
% This creates a master response matrix with all the data entered by the subject and corresponding questions and question types 
% Function was originally written to see how many joins were necessary to obtain all possible data for the responses
% In a real situation, only a part of this query will likely be needed.

% 01/04/07 PJ - flagged incomplete because length of responses has to be
% dynamic to accommodate changing numbers of fields in response table


query_responses =  ['select * from  ' response_table ' ' ...
      'left join `form` on (form.form_id = ' response_table '.form_id) ' ...
      'left join form_x_question on (form_x_question.form_id = form.form_id and ' ...
      'form_x_question.question_iteration = ' response_table '.question_iteration and ' ...
      'form_x_question.form_question_num = ' response_table '.form_question_num) ' ...
      'left join question on (question.question_id = form_x_question.question_id) ' ...
      'left join question_x_data_format on (question_x_data_format.question_id = question.question_id and ' ...
      'question_x_data_format.subquestion = ' response_table '.subquestion) ' ...
      'left join data_format on (data_format.data_format_id = question_x_data_format.answer_format_id) ' ...
      'where ' response_table '.subject_id = ''' subject_id ''' ' ...
      'order by ' response_table '.form_order, form_x_question.form_question_num']
				  
[responses{1:41}] = mysql(query_responses)
			     
return
