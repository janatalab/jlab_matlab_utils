function out_st = ensemble_combine_datastructs(dst,params)
% Combines two Ensemble datastructs in one of the three following ways:
% (1) If the set of variables is identical, the function checks for an 
%     exact match of rows across all variables. If the match is not exact,
%     ensemble_concat_datastruct is called to concatenate the
%     data structures.
% (2) If the two sets of variables differ, a union of the variables will be
%     created if and only if all of the rows of the originally intersecting
%     variable sets match.
%
% (3) If params.heuristic = 'merge', and there is a single intersecting
%     variable, variables from data structure 2 are added to data structure
%     1, populating rows in data structure 1 whose values on the
%     intersecting variable match.
%
% See also: ensemble_concat_datastruct

% 26Jan2013 Petr Janata
% 22Aug2014 PJ - added merge heuristic

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
if isfield(params,'heuristic')
  switch params.heuristic
    case 'merge'
      % Get the intersecting variable
      intersectVars = intersect(dst{1}.vars, dst{2}.vars);
      numIntersect = length(intersectVars);
      if numIntersect > 1
        error('Too many intersecting variables')
      end
      intersectVars = intersectVars{1};
      d1cols = set_var_col_const(dst{1}.vars);
      d2cols = set_var_col_const(dst{2}.vars);
      
      % Merge happens from right to left
      [matchMask, srcIdxs] = ismember(dst{1}.data{d1cols.(intersectVars)}, dst{2}.data{d2cols.(intersectVars)});
      
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
        end
        dst{1}.data{d1cols.(currVar)}(matchMask) = dst{2}.data{d2cols.(currVar)}(srcIdxs(matchMask));
      end
      out_st = dst{1};
  end
end % if isfield (params.'heuristic')

return