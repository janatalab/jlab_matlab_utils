function newStr = reformat_matlab_call(origStr, swap)
% Substitutes context-specific parameter values for parameters specific in
% condition_matlab and stimulus_matlab callbacks in an Ensemble experiment

% 26Jun2014 Petr Janata

% Strip out function parameters
tokens = regexp(origStr,'\''(\w+)\''','tokens');
tokens = [tokens{:}];
pv_pairs = reshape(tokens,2,length(tokens)/2)';

% Deal with variable substitutions as we build up the matlab function
% call
param_str = '';
for iparam = 1:size(pv_pairs,1)
  currParam = pv_pairs{iparam,1};
  switch currParam
    case {'params','return_type'}
      param_str = [param_str sprintf('''%s'',''%s''', currParam, pv_pairs{iparam,2})];
    case 'formid'
      param_str = [param_str sprintf('''%s'', %d', currParam, swap.form_id)];
    case 'subid'
      param_str = [param_str sprintf('''%s'', ''%s''', currParam, swap.subject_id)];
    case 'sessid'
      param_str = [param_str sprintf('''%s'', %d', currParam, swap.session_id)];
    case 'lastvisited'
      param_str = [param_str sprintf('''%s'', %d', currParam, swap.last_visited)];
  end % switch currParam
  
  % Add a comma if needed
  if iparam < size(pv_pairs,1)
    param_str = [param_str ', '];
  end
end % for iparam

% Now reconstruct the MATLAB call
newStr = cell2str([regexp(origStr,'^\w+(','match') param_str ');']);

return