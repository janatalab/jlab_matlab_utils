function out_st = ensemble_combine_datastructs(data_st_array,params)
% Combines two Ensemble datastructs in one of the two following ways:
% (1) If the set of variables is identical, the function checks for an 
%     exact match of rows across all variables. If the match is not exact,
%     ensemble_concat_datastruct is called to concatenate the
%     data structures.
% (2) If the two sets of variables differ, a union of the variables will be
%     created if and only if all of the rows of the originally intersecting
%     variable sets match.
%
% See also: ensemble_concat_datastruct

% 26Jan2013 Petr Janata

if ~iscell(data_st_array) || length(data_st_array) ~= 2
  error('%s: Two data structures are required', mfilename)
end

% Determine whether the sets of variables differ
diffVars = setxor(data_st_array{1}.vars, data_st_array{2}.vars);

if isempty(diffVars) % Sets of variables are the same
  % See if the two structures are the same
  diffRows = setdiff(data_st_array{1}.data, data_st_array{2}.data,'rows');
  
  if ~isempty(diffRows)
    out_st = ensemble_concat_datastruct(data_st_array);
  else
    fprintf('%s: Input data structures are the same. Nothing combined.\n', mfilename);
  end
else % Create a union of the variables if intersection matches
  % Get the set of intersecting variables
  intersectVars = intersect(data_st_array{1}.vars, data_st_array{2}.vars);
  numIntersect = length(intersectVars);
  
  % Get indices into data structure for each data structure
  ds_idxs = cell(1,2);
  for ids = 1:2
    [~,ds_idxs{ids}] = ismember(intersectVars, data_st_array{ids}.vars);
  end
  
  % Compare each of the intersecting variables
  varDiffers = false(1,numIntersect);
  for ivar = 1:numIntersect
    varDiffers(ivar) = ~isempty(setxor(data_st_array{1}.data{ds_idxs{1}(ivar)}, data_st_array{2}.data{ds_idxs{2}(ivar)}));
  end
  
  if any(varDiffers)
    fprintf('%s: Data structures differ in variables and rows. Nothing combined.\n', mfilename);
    out_st = [];
  else
    % Initialize the output data struct with the first data struct
    out_st = data_st_array{1};
    
    % Determine which of the variables from the second data struct need to
    % be copied
    copyVars = setdiff(data_st_array{2}.vars, intersectVars);
    
    % Append the variables that need to be copied to the output data
    % structure list
    out_st.vars = [out_st.vars, copyVars];
    
    % Get the locations of the variables to be copied from the second data
    % struct
    [~,srcCols] = ismember(copyVars, data_st_array{2}.vars);

    % Copy the data struct 2 data to the output structure
    out_st.data = [out_st.data data_st_array{2}.data(srcCols)];
  end
end

return