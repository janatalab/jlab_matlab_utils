function [panas] = extract_panas(expt, sub_list, resp_table, conn_id)
% [panas] = extract_panas(expt, sub_list, resp_table, conn_id)
%
% Extracts PANAS survey data for a given experiment and list of subjects
%
% panas.adj - terms in inventory
% panas.data - nsub X nadj X niter matrix, where niter refers to the number of
% times during the experiment that the survey was administered
%

% 03/04/05 PJ

global SQL_HOST

if nargin < 2
  sub_list = {};
end

try conn_id(1);
catch   
  mysql_make_conn;
  conn_id = 0;
end

% Identify which form_id is for the PANAS Survey
[form_info{1:2}] = mysql(conn_id,'select form_id, form_name from form');
form_id = form_info{1}(strmatch('PANAS Survey',form_info{2},'exact'));


% Get the question IDs that correspond to this form -- in this case a single
% question
quest_id = mysql(conn_id,sprintf('select question_id from form_x_question where form_id=%d',form_id));


% Get the list of PANAS adjectives
[panas.adj, sq_id] = mysql(conn_id,sprintf('select heading, subquestion from question_x_data_format where question_id=%d and subquestion>1', quest_id));


% Pull all of the data for all of the subjects all at once
extract_vars = {'subject_id', 'form_order', 'subquestion', 'response_enum'};
sql_str = sprintf(['select %s from %s where question_id=%d and subquestion>1'], ...
    cell2str(extract_vars,','), resp_table, quest_id);
[panas_resp{1:length(extract_vars)}] = mysql(conn_id,sql_str);

form_col = strmatch('form_order',extract_vars,'exact');
sq_col = strmatch('subquestion',extract_vars,'exact'); 
resp_col = strmatch('response_enum',extract_vars,'exact');

% Determine how many times the PANAS was given, and what the form numbers were
iter_ids = unique(panas_resp{form_col});
niter = length(iter_ids);

% Extract data for each iteration and subject
if isempty(sub_list)
  sub_list = unique(panas_resp{1});
end

nsub_proc = length(sub_list);
panas.subids = sub_list;
panas.data = zeros(nsub_proc,length(sq_id),niter)+NaN;

for iiter = 1:niter
  iter_mask = panas_resp{form_col} == iter_ids(iiter); % iteration mask
  for isub = 1:nsub_proc
    % Make a mask for this subject
    sub_mask = ...
	ismember(panas_resp{strmatch('subject_id',extract_vars,'exact')},sub_list{isub});
    
    % Extract data into final table
    curr_idxs = find(iter_mask&sub_mask);
    panas.data(isub,panas_resp{sq_col}(curr_idxs)-1,iiter) = log2(panas_resp{resp_col}(curr_idxs));
    
  end % for isub
end % for iiter

if ~conn_id
  mysql(conn_id,'close')
end
return
