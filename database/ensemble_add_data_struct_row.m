function indata = ensemble_add_data_struct_row(indata,varargin)

% adds a given row of data to a given ensemble data struct
% 
% REQUIRES
%   indata - ensemble data struct
%   varargin - key/value pairs representing a row of data
%       'key' - always a string matching a column name in indata.vars
%       'value' - could be any type, is added as the value of 'key'
% 
% NOTE: all values are appended to the end of each existing row
% NOTE: fills out short columns with empty values of the appropriate type
% 
% FIXME: allow addition of multiple rows? pass in vectors instead of scalar
% or single-item types? how would we differentiate vectors for rows and
% vectors that should be an individual indexed item in indata? any varargin
% switch or keyword would have to be something that would never ever ever
% be a variable name in an ensemble data struct...
% 
% FB 2009.03.09
% PJ 2010.11.08 - fixed handling of single items in that were not cell arrays and
%                 which generated an excess number of rows by virtue of
%                 their length.

if ~is_ensemble_datastruct(indata)
  error('indata is not a valid ensemble data structure\n');
end

% get var names
cols = set_var_col_const(indata.vars);

% add row
for iv = 1:2:nargin-1
  fld = varargin{iv};
  if isfield(cols,fld)
%     if iscell(indata.data{cols.(fld)})
%       indata.data{cols.(fld)}{end+1} = varargin{iv+1};
%     else
%       inclass = class(indata.data{cols.(fld)});
%       vclass  = class(varargin{iv+1});
%       if isempty(strmatch(inclass,vclass))
%         error(['data struct column type (%s) and value type (%s) for '...
%             '%s key do not match'],inclass,vclass,fld);
%       else
        indata.data{cols.(fld)} = [indata.data{cols.(fld)}; varargin{iv+1}];
%         indata.data{cols.(fld)}(end+1) = varargin{iv+1};
%         indata.data{cols.(fld)} = [
%       end
%     end
  end % if isfield(cols,
end % for ivar =

% find maximum column length
%maxn = max(cellfun(@length,indata.data));
dims = cellfun(@size,indata.data,'UniformOutput',false);
dims = cat(1,dims{:});
maxn = max(dims(:,1));

% fill out short columns
for id = 1:length(indata.data)
  while length(indata.data{id}) < maxn
    if iscell(indata.data{id})
      if isempty(indata.data{id})
        indata.data{id} = cell(maxn,1);
      else
        indata.data{id} = [indata.data{id}; ...
            cell(maxn-length(indata.data{id}),1)];
      end
    else
      if isempty(indata.data{id})
        indata.data{id} = nan(maxn,1);
      else
        indata.data{id} = [indata.data{id}; ...
            nan(maxn-length(indata.data{id}),1)];
      end
    end
  end
end
