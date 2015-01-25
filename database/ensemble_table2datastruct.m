function data_st = ensemble_table2datastruct(data_t)
% Converts a MATLAB table variable into an Ensemble datastruct

% 26Jan2015 Petr Janata

data_st = ensemble_init_data_struct;
data_st.vars = data_t.Properties.VariableNames;
nvars = length(data_st.vars);

data_st.data = cell(1,nvars);
data_st.name = data_t.Properties.Description;

for ivar = 1:nvars
    data_st.data{ivar} = data_t{:,ivar};  
end

end