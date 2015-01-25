function out_st = ensemble_reshape_data(data_st,params)
% Reshapes a data struct by making new variables out of the unique levels
% in xfmVar.
%
% out_st = ensemble_reshape_data(data_st,params);
% 
% USAGE:
%   data_st - an Ensemble data struct
%   params - see below
%
% In order for this tranformation to work it requires 3 sets of variables,
% specified in params. The following fields must be present in
% params.ensemble_reshape_data_st
%
% REQUIRED:
% xfmVar - This is the name of the variable whose unqiue values will form
%          new variables in the out_st Ensemble data struct. This will
%          typically be a 'compqid_str' variable (created by
%          ensemble_attach_compqid_str)
%
% valueVars - this is the variable whose values will populate the new
%          variables. This would typically be 'response_enum' or
%          'response_text'. If multiple values are specified, both are
%          checked. If a value is contained in both of them, an error is
%          thrown. If a value is contained in only 1, it is copied. If
%          neither contains a value, a NaN or empty string is placed,
%          depending on type.
%
% keyVars - this is the set of variables, which in association with xfmVar,
%          will identify single (unique) rows in the data struct. If unique
%          rows are not found, and error will be generated
%
% OPTIONAL:
% var_name_map - If present, this structure is used to map levels of xfmVar 
%          to names that will be used as variable names in the output
%          data_st
% copyVars - additional variables from the data_st that will be copied to
%          the output structure
%
% See also: ensemble_data_by_question(), ensemble_export_respnstim(), ensemble_vals2vars()
%
% The difference between ensemble_reshape_data_st and ensemble_vals2vars is
% that the latter will not handle multi-variable keys.

% 04Sep2014 Petr Janata

% Turn off categorical class warnings
warning('off','stats:categorical:subsasgn:NewLevelsAdded')

%% Check input parameters
if nargin < 2
  error('%s: data_st and params inputs required', mfilename)
end

if ~isfield(params,mfilename)
  error('params.%s structure is required',mfilename)
end

requiredVars = {'xfmVar','valueVars','keyVars'};
missingMask = ~ismember(requiredVars,fieldnames(params.(mfilename)));
if any(missingMask)
  error('Required fields missing from params.%s: %s', mfilename, cell2str(requiredVars(missingMask),','))
end

if isfield(params,'verbose')
  verbose = params.verbose;
else
  verbose = 1;
end

% Make variable lists a bit more accessible
xfmVar = params.(mfilename).xfmVar;

keyVars = [params.(mfilename).keyVars xfmVar];
nkeys = length(keyVars);

valueVars = params.(mfilename).valueVars;

if isfield(params.(mfilename),'copyVars')
  copyVars = params.(mfilename).copyVars;
else
  copyVars = {};
end
ncopy = length(copyVars);

% Get column indexing for input data
cols = set_var_col_const(data_st.vars);

%% Perform any desired filtering
if isfield(params,'filt')
  data_st = ensemble_filter(data_st, params.filt);
end

%% Extract the key variables and xfmVar and check for uniqueness
nrows = size(data_st.data{cols.(keyVars{1})},1);

% Convert each variable to a categorical type
for ikey = 1:nkeys
  keymtx(:,ikey) = nominal(data_st.data{cols.(keyVars{ikey})});
end

% Check for uniqueness
nonUniqueMask = check_unique_rows(keymtx, verbose);

if any(nonUniqueMask)
  error('Key variables do not return unique values')
end

%% Get the set of output keys (all keys but xfmVar key)
outkeymtx = unique(keymtx(:,1:end-1),'rows');
noutRows = size(outkeymtx,1);

%% Extract a matrix with all possible values
valmtx = data_st.data(ismember(data_st.vars,valueVars));
nvals = length(valmtx);

% Get types for each of the value variables
valtype = cell(1,nvals);
for ival = 1:nvals
  valtype{ival} = class(valmtx{ival});
end

%% Initialize the output structure
out_st = ensemble_init_data_struct;

% Add the variables we want to copy to the list and initialize their types
if ncopy
  out_st.vars = copyVars;
  ocols = set_var_col_const(out_st.vars);
  for icopy = 1:ncopy
    currVar = copyVars{icopy};
    currType = class(data_st.data{cols.(currVar)});
    switch currType
      case {'numeric','double'}
        out_st.data{ocols.(currVar)} = nan(noutRows,1);
      case 'logical'
        out_st.data{ocols.(currVar)} = false(noutRows,1);
      case 'cell'
        out_st.data{ocols.(currVar)} = cell(noutRows,1);
      otherwise
        error('No initialization for type: %s', currType)
    end
  end
end

%% Create new variables in output structure
[levelMaskMtx, newVars] = make_mask_mtx(data_st.data{cols.(xfmVar)});
numNew = length(newVars);
newSrc = cell(1,numNew);

for inew = 1:numNew
  currLevel = newVars{inew};
  
  % Find rows in the input data corresponding to this variable
  levelMask = levelMaskMtx(:,strcmp(newVars,currLevel));

  % See if we want to remap its name
  if isfield(params.(mfilename).var_name_map,currLevel)
    varName = params.(mfilename).var_name_map.(currLevel);
  else
    varName = currLevel;
  end
  
  % Place the name in the output structure if necesary
  out_st.vars{end+1} = varName;
  ocols = set_var_col_const(out_st.vars);
 
  % Check which of the value variables we have data in for this level
  haveData = false(1,nvals);
  for ival = 1:nvals
    currData = data_st.data{cols.(valueVars{ival})}(levelMask);
    switch valtype{ival}
      case {'numeric','double','logical'}
        haveData(ival) = any(currData);
      case 'cell'
        haveData(ival) = any(~cellfun('isempty', currData));
    end
  end
  
  if ~any(haveData)
    error('No data available in any of the value variables')
  end
  
  if sum(haveData) > 1
    error('More than one value variable has data for level (%s): %s', ...
      currLevel, cell2str(valueVars(haveData),','))
  end
  
  valVar = valueVars{haveData};
  currType = valtype{haveData};
  newSrc{inew} = valVar;
  
  % Initialize the output variable
  switch currType
    case {'numeric','double'}
      out_st.data{ocols.(varName)} = nan(noutRows,1);
    case 'logical'
      out_st.data{ocols.(varName)} = false(noutRows,1);
    case 'cell'
      out_st.data{ocols.(varName)} = cell(noutRows,1);
    otherwise
      error('No initialization for type: %s', currType)
  end
end % new variable initialization


%% Copy data to the new variable

% Loop over values in the key matrix
for irow = 1:noutRows
  currKey = outkeymtx(irow,:);
  outkeyMask = ismember(keymtx(:,1:end-1),currKey,'rows');
  
  for inew = 1:numNew
    currLevel = newVars{inew};
    
    % See if we want to remap its name
    if isfield(params.(mfilename).var_name_map,currLevel)
      varName = params.(mfilename).var_name_map.(currLevel);
    else
      varName = currLevel;
    end
 
    % Find rows corresponding to this variable
    levelMask = levelMaskMtx(:,strcmp(newVars,currLevel));
    
    % Get the composite mask
    compMask = outkeyMask & levelMask;
    
    % Skip row if not relevant
    if ~any(compMask)
      continue
    end
    
    % Copy the value variable
    out_st.data{ocols.(varName)}(irow) = data_st.data{cols.(newSrc{inew})}(compMask);
  
 end % for inew
 
 % Now copy the other variables we want to carryover
 for icopy = 1:ncopy
   tmp = data_st.data{cols.(copyVars{icopy})}(outkeyMask);
   if ~iscell(tmp) && any(isnan(tmp))
     continue
   end
   
   outval = unique(tmp);
   if numel(outval) > 1
     error('Too many (%d) values found for variable we want to copy (%s)', numel(outval), copyVars{icopy})
   else
     out_st.data{ocols.(copyVars{icopy})}(irow) = outval;
   end
 end % for icopy
    
 
end % for irow


end % function