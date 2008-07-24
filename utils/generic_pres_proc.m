function [data,meta] = generic_pres_proc(fname,p)
% [data,meta] = generic_pres_proc(fname,params)
%
% fname - Presentation file to parse
% params - a structure that controls details of the parse and what meta data
%          should be returned
%
% Options for params structure
% 
% 'skipline','delim','verbose' - control loadtxt() behavior
% 'strip_colnames' - if true, the first row of data containing column names is stripped
% 'mask' - a structure array specifying if and how any masks should be created
%
% Fields within the  mask structure
% 'name' - name of the mask which becomes the field name in the output structure
% 'crit_col_lbl' - name in the PL (column indexing) structure to be used
% 'crits' - criteria (events) to match
% 'is_partial' - is the value in crits a partial search string

% 12/25/06 Petr Janata
data = [];
meta = [];

if nargin < 2
  p = struct([]);
end

p = proc_input_params(p);
  
if ~exist(fname)
  fprintf('Did not find file: %s\n', fname);
  return
end

% Load the data from the Presentation file using the EEGLAB function loadtxt
data = loadtxt(fname, 'skipline', p.skipline, 'delim', p.delim, 'verbose', p.verbose);
col_names = data(1,:);

% Map the column names to column index values
PL = set_pres_col_const(col_names);
meta.PL = PL;

% Strip the row of column names
if p.strip_colnames
  data(1,:) = [];
end

% Create masks if desired
nmasks = length(p.mask);
for imask = 1:nmasks
  mask_name = p.mask(imask).name;
  crit_col = p.mask(imask).crit_col_lbl;
  mask_crit = p.mask(imask).crits;
  if ~iscell(mask_crit)
    mask_crit = {mask_crit};
  end
  
  tmpdata = data(:,PL.(crit_col));
  % We have to set up some handling for mixed string and numeric types
  str_mask = cellfun(@isstr,tmpdata);
  num_mask = cellfun(@isnumeric,tmpdata);
  
  is_crit_str = cellfun(@isstr, mask_crit);
  is_crit_num = cellfun(@isnumeric, mask_crit);
  
  if ~(all(is_crit_str) | all(is_crit_num))
    fprintf(['generic_pres_proc: Cannot handle mixed type criteria right now' ...
	  ' ...\n']);
    continue
  end
  
  meta.(mask_name) = false(size(tmpdata,1),1);
  if all(is_crit_str)
    if ~p.mask(imask).is_partial;
      meta.(mask_name)(str_mask) = ismember(tmpdata(str_mask), mask_crit);
    else
      items = strfind(tmpdata(str_mask), mask_crit{1});
      meta.(mask_name)(str_mask) = ~cellfun(@isempty,items);
    end
  else
    meta.(mask_name)(num_mask) = ismember([tmpdata{num_mask}], [mask_crit{:}]);
  end
  
end % for imask
  
end

function p = proc_input_params(p)
  % Parameters for loadtxt
  if ~isfield(p,'skipline') p.skipline = 3; end
  if ~isfield(p,'delim') p.delim = 9; end
  if ~isfield(p,'verbose') p.verbose = 'off'; end

  % Post-processing of loaded file
  if ~isfield(p,'strip_colnames') p.strip_colnames = true; end
  
  % Mask creation
  if ~isfield(p,'mask') p.mask = []; end
  
end
