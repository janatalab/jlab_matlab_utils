function ds = ensemble_datastruct2dataset(data_st)
% Converts an Ensemble data struct to a dataset variable

% 14Sep2014 Petr Janata

ds = dataset(data_st.data{:},'VarNames',data_st.vars);

end