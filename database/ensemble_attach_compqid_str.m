function data_st = ensemble_attach_compqid_str(data_st,~)
% Attaches a compqid_str variable to the Ensemble data struct if none
% exists.
%
% INPUT:
%   data_st - must contain at least question_id and subquestion variables
%
% OUTPUT:
%   Items in the compqid_str variable have the form
%   s<question_id>_<subquestion>, e.g. s1032_01

% 04Sep2014 Petr Janata

if ~all(ismember({'question_id','subquestion'},data_st.vars))
  error('%s: Both question_id and subquestion variables required', mfilename)
end

cols = set_var_col_const(data_st.vars);

% Create the initial string
tmp = strcat('s', num2str(data_st.data{cols.question_id}),'_',num2str(data_st.data{cols.subquestion},'%02d'));

% Place into a cell array of strings
tmp = cellstr(tmp);

% Remove whitespace
tmp = regexprep(tmp,'\s','');

data_st.vars{end+1} = 'compqid_str';
data_st.data{end+1} = tmp;

end