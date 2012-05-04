function data_st = ensemble_check_compqid(data_st)
% Checks to see if the data structure has a composite question ID variable.
%
% If it has none, it checks to see if there are question_id and subquestion
% variables from which to construct a compqid variable, and proceeds to do so if
% these are available. Otherwise, the function returns empty.
%
% data_st = ensemble_check_compqid(data_st);
%

% 02/06/07 Petr Janata

col = set_var_col_const(data_st.vars); % set the column constants

if isfield(col,'compqid')
  return % nothing to construct, so return what we received
end

if ~isfield(col,'compqid') && ~all(isfield(col,{'question_id','subquestion'}))
  fprintf(['%s: Did not find necessary question and subquestion or composite question ID' ...
	' information in the input data'], mfilename);
  data_st = [];
  return
end

% Build the compqid variable if necessary
data_st.vars{end+1} = 'compqid';
data_st.data{end+1} = make_compqid(data_st.data{col.question_id}, ...
    data_st.data{col.subquestion});
  
end
