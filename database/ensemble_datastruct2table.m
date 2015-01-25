function ds = ensemble_datastruct2table(data_st)
% Converts an Ensemble data struct to a table variable

% 17Sep2014 Petr Janata

ds = table(data_st.data{:},'VariableNames',data_st.vars);
ds.Properties.Description = data_st.name;

end