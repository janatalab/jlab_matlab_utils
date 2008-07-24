function FD = set_form_col_const(var_list);
% Figures out which data columns correspond to different database fields in the form table.
%
% Use of the FD structure when indexing into form data ensures portability of code
% in the event that the field ordering in the database changes or if only a subset
% of fields are retrieved in a query.
%
% var_list is a cell string containing field names
% FD is a structure with constants (column ids) assigned for the following
% fields:
%
% SUB_ID -- subject ID
% STIM_ID -- stimulus ID
% DATE_TIME -- datestamp
% QUEST_ID -- question ID
% QUEST_TXT -- question text
% SUBQUEST_TXT -- subquestion text (heading)
% SUBQUEST_ID -- subquestion ID
% QUEST_ITER -- question iteration
% RESP_ENUM -- response value
% RESP_TXT -- response text
% RESP_ID -- incrementing response ID counter
% SESS_ID -- session ID
% FORM_ID -- form ID
% 

% March 2005, Petr Janata
% 06/22/06 PJ - added session_id handling
% 01/04/07 PJ - added form_id handling

nvars = length(var_list);

FD = struct('SUB_ID',[], ...
    'STIM_ID',[], ...
    'DATE_TIME',[], ...
    'QUEST_ID',[], ...
    'QUEST_TXT',[], ...
    'SUBQUEST_TXT',[], ...
    'SUBQUEST_ID',[], ...
    'QUEST_ITER',[], ...
    'RESP_ENUM',[], ...
    'RESP_TXT',[], ...
    'RESP_ID',[], ...
    'SESS_ID',[], ...
    'FORM_ID',[] ...
    );

for ivar = 1:nvars
  switch var_list{ivar}
    case 'subject_id'
      FD.SUB_ID = ivar;
      
    case 'stimulus_id'
      FD.STIM_ID = ivar;
      
    case 'date_time'
      FD.DATE_TIME = ivar;
      
    case 'question_id'
      FD.QUEST_ID = ivar;
      
    case 'subquestion'
      FD.SUBQUEST_ID = ivar;

    case 'question_iteration'
      FD.QUEST_ITER = ivar;
    
    case 'response_enum'
      FD.RESP_ENUM = ivar;
    
    case 'response_text'
      FD.RESP_TXT = ivar;
    
    case 'question_text'
      FD.QUEST_TXT = ivar;
    
    case 'heading'
      FD.SUBQUEST_TXT = ivar;
      
    case 'response_id'
      FD.RESP_ID = ivar;
  
    case 'session_id'
      FD.SESS_ID = ivar;
    
    case 'form_id'
      FD.FORM_ID = ivar;
  end % switch
end % for ivar

return
