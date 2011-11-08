function data_st = ensemble_remove_vars_from_datastruct(data_st, vars2remove)
% Removes variables specified in vars2remove cell array of strings from the
% Ensemble datastruct data_st

% 08Nov2011 Petr Janata

rmMask = ismember(data_st.vars, vars2remove);

data_st.vars(rmMask) = [];
data_st.data(rmMask) = [];

return