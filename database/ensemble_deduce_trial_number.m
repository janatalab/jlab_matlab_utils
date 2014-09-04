function out_st = ensemble_deduce_trial_number(resp_st,params)
% Deduces trial numbers associated with repeated forms in a response table
% structure. Returns trial numbers (starting at 1 for each subject
% indicated in the subject_id field) in a new field called trial_number. 
%
% out_st = ensemble_deduce_trial_number(resp_st,params);
%
% This function is useful for generating a variable that can be used to
% group multiple rows that belong together conceptually. Having such
% information is necessary for rearranging functions that convert between
% long and wide formats, such as when different question_id or compqid
% values are broken out into their own variables. For example,
% export_respnstim accomplishes this rearrangement by virtue of a
% consistent stimulus_id across a set of rows. However, non-stimulus
% related blocks of questions have no grouping information, and so
% trial_number fulfills that purpose.
%
% 

% 03Sep2014 Petr Janata

% Make sure we have the columns we need
requiredVars = {'subject_id','response_order','form_order'};
if ~all(ismember(requiredVars, resp_st.vars))
  error('Required variables (%s) not in data struct', cell2str(requiredVars,','))
end

% Add the new column
resp_st.vars{end+1} = 'trial_number';
resp_st.data{end+1} = zeros(size(resp_st.data{1},1),1);

% Get column names and indices
rcols = set_var_col_const(resp_st.vars);

% Get the number of subjects
subids = unique(resp_st.data{rcols.subject_id});
nsubs = length(subids);

% Generate a matrix of response and form orders
rofo = [resp_st.data{rcols.response_order} resp_st.data{rcols.form_order}];

for isub = 1:nsubs
  currSub = subids{isub};
  submask = strcmp(resp_st.data{rcols.subject_id}, currSub);
  
  newro = [1; diff(rofo(submask,1))>0];
  newsamefo = [1; diff(rofo(submask,2))<1];
  
  rofomask = newro & newsamefo;
  
  resp_st.data{rcols.trial_number}(submask) = cumsum(rofomask);
  
end % for isub

return