function out_st = ensemble_combine_datastructs(dst,params)
% Combines two Ensemble datastructs in one of the three following ways:
% (1) If the set of variables is identical, the function checks for an 
%     exact match of rows across all variables. If the match is not exact,
%     ensemble_concat_datastruct is called to concatenate the
%     data structures.
%
% (2) If the two sets of variables differ, a union of the variables will be
%     created if and only if all of the rows of the originally intersecting
%     variable sets match.
%
% (3) If params.heuristic = 'merge', variables from data structure 2 are
%     added to data structure 1, populating rows in data structure 1 whose
%     values on the set of intersecting variables match. For each
%     structure, uniqueness of rows (created from intersecting variables)
%     is required and tested.
%
% In order to accomplish the multiple variable merge, it was necessary to
% utilize the categorial data types from the Statistics Toolbox.
%
% See also: ensemble_concat_datastruct

% 26Jan2013 Petr Janata
% 22Aug2014 PJ - added merge heuristic
% 03Sep2014 PJ - added support for merge based on multiple intersecting
%                variables

if ~iscell(dst) || length(dst) ~= 2
  error('%s: Two data structures are required', mfilename)
end

% Determine whether the sets of variables differ
diffVars = setxor(dst{1}.vars, dst{2}.vars);

if isempty(diffVars) % Sets of variables are the same
  % See if the two structures are the same
  diffRows = setdiff(dst{1}.data, dst{2}.data,'rows');
  
  if ~isempty(diffRows)
    out_st = ensemble_concat_datastruct(dst);
  else
    fprintf('%s: Input data structures are the same. Nothing combined.\n', mfilename);
  end
else % Create a union of the variables if intersection matches
  % Get the set of intersecting variables
  intersectVars = intersect(dst{1}.vars, dst{2}.vars);
  numIntersect = length(intersectVars);
  
  % Get indices into data structure for each data structure
  ds_idxs = cell(1,2);
  for ids = 1:2
    [~,ds_idxs{ids}] = ismember(intersectVars, dst{ids}.vars);
  end
  
  % Compare each of the intersecting variables
  varDiffers = false(1,numIntersect);
  for ivar = 1:numIntersect
    % Make sure that if we are dealing with a cell array that contains
    % strings that any empty values are empty strings, rather than empty
    for ist = 1:2
      if iscell(dst{ist}.data{ds_idxs{ist}(ivar)})
        emptyMask = cellfun('isempty',dst{ist}.data{ds_idxs{ist}(ivar)});
        [dst{ist}.data{ds_idxs{ist}(ivar)}{emptyMask}] = deal('');
      end
    end
    varDiffers(ivar) = ~isempty(setxor(dst{1}.data{ds_idxs{1}(ivar)}, dst{2}.data{ds_idxs{2}(ivar)}));
  end
  
  if any(varDiffers)
    if ~isfield(params,'heuristic')
      fprintf('%s: Data structures differ in variables and rows. Nothing combined.\n', mfilename);
    end
    out_st = [];
  else
    % Initialize the output data struct with the first data struct
    out_st = dst{1};
    
    % Determine which of the variables from the second data struct need to
    % be copied
    copyVars = setdiff(dst{2}.vars, intersectVars);
    
    % Append the variables that need to be copied to the output data
    % structure list
    out_st.vars = [out_st.vars, copyVars];
    
    % Get the locations of the variables to be copied from the second data
    % struct
    [~,srcCols] = ismember(copyVars, dst{2}.vars);

    % Copy the data struct 2 data to the output structure
    out_st.data = [out_st.data dst{2}.data(srcCols)];
  end
end

% See if the out_st is still empty. If so, check whether a different
% heurisitic has been specified.
CAN_HANDLE_MULTIPLE_VARS = 1;
if isfield(params,'heuristic')
  switch params.heuristic
    case 'merge'
      % Get the intersecting variable
      intersectVars = intersect(dst{1}.vars, dst{2}.vars);
      
      numIntersect = length(intersectVars);
      if numIntersect > 1 && ~CAN_HANDLE_MULTIPLE_VARS
        error('Too many intersecting variables')
      end
      
      % Create tables for each data structure consisting only of the
      % intersecting variables. We can the test these tables for
      % uniqueness of rows, and find row matches between the two data
      % structs.
      % To accomplish this, we have to first convert the data to a
      % categorical class (requires Statistics Toolbox), which each of the
      % variables being a nominal type. Unique rows can then be found.
      ordinalTbl = cell(1,2);
      for ids = 1:2
        [~,idxs] = ismember(intersectVars,dst{ids}.vars);
        tbl{ids} = dst{ids}.data(idxs);
        for iint = 1:length(intersectVars)
          ordinalTbl{ids}(:,iint) = ordinal(tbl{ids}{iint});
        end
        
        % Check for uniqueness of rows
        uniqueRows = unique(ordinalTbl{ids},'rows');
        numNonUnique = length(ordinalTbl{ids}) - length(uniqueRows);
        if numNonUnique
          fprintf('Found %d non-unique rows in data struct: %s\n',numNonUnique,dst{1}.name);
          % Find the non-unique rows
          [~,idxs] = ismember(ordinalTbl{ids},uniqueRows,'rows');
          
          % Tabulate idxs to see which appears more than once
          t = tabulate(idxs);
          
          % Display the nonUnique rows
          uniqueRows(t(:,2)>1,:)
          
          fprintf('Cannot merge data ...\n');
          return
        end       
      end
      
      % Convert our ordinal table to a nominal table in order to accomplish
      % the following ismember operation
      [matchMask, srcIdxs] = ismember(nominal(ordinalTbl{1}),nominal(ordinalTbl{2}),'rows');
      
      d1cols = set_var_col_const(dst{1}.vars);
      d2cols = set_var_col_const(dst{2}.vars);
      
      % Figure out which variables we'll be appending to dst1
      copyVars = setdiff(dst{2}.vars, dst{1}.vars);
      dst{1}.vars = [dst{1}.vars copyVars];
      d1cols = set_var_col_const(dst{1}.vars);
      nvars = length(copyVars);
      for ivar = 1:nvars
        currVar = copyVars{ivar};
        
        % Initialize the d1 data cell to the correct type
        if iscell(dst{2}.data{d2cols.(currVar)}(1))
          dst{1}.data{d1cols.(currVar)} = cell(size(matchMask));
        elseif isnumeric(dst{2}.data{d2cols.(currVar)}(1))
          dst{1}.data{d1cols.(currVar)} = nan(size(matchMask));
        elseif islogical(dst{2}.data{d2cols.(currVar)}(1))
          dst{1}.data{d1cols.(currVar)} = false(size(matchMask));
        end
        dst{1}.data{d1cols.(currVar)}(matchMask) = dst{2}.data{d2cols.(currVar)}(srcIdxs(matchMask));
      end
      out_st = dst{1};
    otherwise
      fprintf('%s: Unknown heuristic: %s\nData structures not combined ...\n', mfilename, params.heuristic);
  end
end % if isfield (params.'heuristic')

return