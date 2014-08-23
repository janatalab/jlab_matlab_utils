function [out_st, newVars] = ensemble_vals2vars(data_st, params)
% Rearranges an ensemble data struct such that unique combinations of
% values in one or more variables specified in params.src_vars become their
% own variables in the output structure.
%
% Additional parameters that govern how this function work are:
% params.by_var - Unique values of this variable are looped over when
%                 creating the rearranged matrix. This variable should have
%                 exactly one instance of each unique combination of the
%                 values of the source variables, i.e. new variable. If
%                 this field is not specified, the first variable in the
%                 data structure variables is used. 
%
% params.value_var - name of the variable that contains the actual values
%                    that will be placed in the new variables. If this is
%                    not specified, 'value' is used. 
% params.carryover_vars - list of variable names from original data
%                  structure to carry over to the new structure.

% 09Nov2011 Petr Janata
% 22Aug2014 PJ - fixed handling of non-cell by_var; added carryover_vars
%                option

if ~isfield(params,'src_vars')
	error('%s: Variables whose values to transform not specified in params.src_vars', mfilename)
end

xfmVars = params.src_vars;
if ~iscell(xfmVars)
  xfmVars = {xfmVars};
end
nvars = length(xfmVars);

srcVars = data_st.vars;
srcCols = set_var_col_const(data_st.vars);
[varMask,srcIdxs] = ismember(xfmVars, srcVars);
if ~all(varMask)
	error('%s: Desired variables are not all present', mfilename);
end

% Figure out which variable we are going to use for looping when we
% rearrange the data
if isfield(params,'by_var')
	byVar = params.by_var;
else
	byVar = srcVars{1};
end
fprintf('Will loop using variable: %s\n', byVar);

% Figure out which variable we are pulling data from
if isfield(params,'value_var')
	valVar = params.value_var;
else
	valVar = 'value';
end
if ~ismember(valVar, srcVars)
	error('Variable from which values will be rearranged (%s) could not be found', valVar)
end

% Initialize the output data structure
out_st = data_st;
out_st.vars = {};
out_st.data = {};

% Extract the desired variables
vardata = data_st.data(srcIdxs);

% Convert any numeric data to strings
for ivar = 1:nvars
	% If it is already in numeric format, but not in a cell, place in cell so
	% that it can be converted properly to a string without whitespace
	if isnumeric(vardata{ivar}) 
		vardata{ivar} = num2cell(vardata{ivar});
	end
	if ~ischar(vardata{ivar}{1})
		fprintf('Converting numeric values to cell array of strings for variable %s\n', xfmVars{ivar});
		vardata{ivar} = cellfun(@num2str,vardata{ivar},'UniformOutput',false);
	end
end

% Get a list of unique values along each dimension
unique_vals = cell(nvars,1);
for ivar = 1:nvars
	unique_vals{ivar} = unique(vardata{ivar});
end

% Reformat vardata
vardata = cat(2,vardata{:});

% Get the unique combinations that are going to form the new variables
fprintf('Finding unique combinations of values of source variables\n'); 
[newVarMask, uniqueCombos] = make_mask_mtx(vardata);
nNewVars = size(uniqueCombos,1);
newVars = cell(nNewVars,1);
for inew = 1:nNewVars
	newVars{inew,1} = cell2str(uniqueCombos(inew,:),'_');
end

% Sanitize new variable names
newVars = strrep(newVars,'.','p');

%
% Set the new variable list for the output structure
%

% Figure out which variables we are carrying over
if isfield(params,'carryover_vars')
  carryoverVars = params.carryover_vars;
else
  carryoverVars = setdiff(srcVars, xfmVars);
end
nCarryover = length(carryoverVars);

% Set the variable list
out_st.vars = [carryoverVars newVars'];
outcols = set_var_col_const(out_st.vars);
out_st.data = cell(1,length(out_st.vars));

% Loop using unique values of the looping variable
uniqueVals = unique(data_st.data{srcCols.(byVar)});
nvals = length(uniqueVals);

fprintf('Looping over %d values for variable %s\n', nvals, byVar);
for ival = 1:nvals
	if mod(ival,20) == 0
		fprintf('%d', ival);
	else
		fprintf('.');
  end
  if iscell(uniqueVals)
    currVal = uniqueVals{ival};
  else
    currVal = uniqueVals(ival);
  end
	valMask = ismember(data_st.data{srcCols.(byVar)}, currVal);

	% Copy carryover variables
	for icarry = 1:nCarryover
		currVar = carryoverVars{icarry};
		
		if strcmp(currVar,valVar)
			continue
		end		
		
		% Make sure we have non-empty values
		if iscell(data_st.data{srcCols.(currVar)}(1))
			if all(cellfun('isempty',data_st.data{srcCols.(currVar)}(valMask)))
				continue
			end
		end
			
		uniqueOldVal = unique(data_st.data{srcCols.(currVar)}(valMask));
		if length(uniqueOldVal) > 1
			% Check to see if we are dealing with only NANs. If so, simply skip
			% over this variable
			if ~iscell(uniqueOldVal) && all(isnan(uniqueOldVal))
				continue
			else
				error('\nMore than 1 unique old value')
			end
		end
		if iscell(uniqueOldVal)
			uniqueOldVal = uniqueOldVal{1};
    end
    if isnumeric(uniqueOldVal)
      out_st.data{outcols.(currVar)}(ival,1) = uniqueOldVal;
    else
      out_st.data{outcols.(currVar)}{ival,1} = uniqueOldVal;
    end
  end % end of carryover variables
	
	% Now, loop over all the new variables and pull the relevant data
	nNewVars = length(newVars);
	for inew = 1:nNewVars
		% Find matching rows
		currMask = newVarMask(:,inew);
	
		% Get the intersection of the current variable mask and stimulus mask
		compositeMask = valMask & currMask;
		if sum(compositeMask) > 1
			error('\nToo many values for instance of variable %s and variable %s', currVal, newVars{inew})
		end
		
		% Assign the value
		tmpval = data_st.data{srcCols.(valVar)}(compositeMask);
    if isempty(tmpval)
      if isnumeric(out_st.data{outcols.(newVars{inew})}(1))
        tmpval = NaN;
      else
        tmpval = '';
      end
    end
    
		if iscell(tmpval)
			tmpval = tmpval{1};
    end
    if isnumeric(tmpval)
      out_st.data{outcols.(newVars{inew})}(ival,1) = tmpval;
    else
      out_st.data{outcols.(newVars{inew})}{ival,1} = tmpval;
    end
		
	end
end % for ival

fprintf('\n')

% Remove the value variable
out_st = ensemble_remove_vars_from_datastruct(out_st, {valVar});

return
