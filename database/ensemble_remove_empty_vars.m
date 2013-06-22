function data_st = ensemble_remove_empty_vars(data_st)
% Removes variables from an Ensemble datastruct whose values are all empty
%

% 21Jun2013 Petr Janata

nvars = length(data_st.vars);
rmvars = [];
for ivar = 1:nvars
  if iscell(data_st.data{ivar})
    if all(cellfun('isempty',data_st.data{ivar}))
      rmvars(end+1) = ivar;
    end
  else
    if all(isnan(data_st.data{ivar}))
      rmvars(end+1) = ivar;
    end
  end
end

fprintf('Removing %d empty variables from data struct\n', length(rmvars));
data_st.vars(rmvars) = [];
data_st.data(rmvars) = [];

end

