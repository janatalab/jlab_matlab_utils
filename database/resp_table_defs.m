% Database response table column mappings
%
% resp_table_defs.m
%
% The reason for packing them into an RT struct is that some of the same column
% names are used in different tables, but might appear in a different order.

% 03/04/05 PJ

RT.RESP_ID = 1;  % response_id
RT.DATATIME = 2; % date_time
RT.RESP_TXT = 3; % response_text
RT.RESP_ENUM = 4; % response_enum
RT.QUEST_ID = 5; % question_id
RT.FORM_QUEST_NUM = 6; % form_question_num
RT.QUEST_ITER = 7; % question_iteration
RT.SUBQUEST = 8; % subquestion
RT.STIM_ID = 9; % stimulus_id
RT.FORM_ID = 10; % form_id
RT.FORM_ORDER = 11; % form_order
RT.SESS_ID = 12; % session_id
RT.SUBJ_ID = 13; % subject_id
RT.EXPER_ID = 14; % experiment_id
