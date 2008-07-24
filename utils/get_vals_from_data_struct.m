function columnData = get_vals_from_data_struct(data_st,columnName)
%
% [data_st] = ensemble_get_vals_from_data_struct(data_st,column)
%
%
% returns a column of a data struct initialized by
% ensemble_init_data_struct
%
% 14 Feb, 2007 First Version, S.T.

colIdx = strmatch(columnName,data_st.vars);
columnData = data_st.data{colIdx};




