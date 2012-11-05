function st = mkstruct(fields,args)
% Returns a structure using the provided fieldnames and values
%
% st = mkstruct(fields,args);
%
%
%
% INPUTS
% fields - cell array of fieldnames to populate. 
% args - a list of param/value pairs containing the name of the field and the
%        value it should be populated with
%
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author:
% 12/04/06 Petr Janata 

st = [];

if nargin < 2
  args = {};
end

nflds = length(fields);
% Initialize a structure with the fields
for ifld = 1:nflds
  st.(fields{ifld}) = [];
end

% Handle the input arguments
nargs = length(args);
if nargs
  for iarg = 1:2:nargs
    tag = args{iarg};
    if ~isstr(tag)
      errorstr = sprintf('Parameter tag must be a string');
      error(errorstr);
    end
    
    % Make sure the tag exists
    if ~any(ismember(lower(fields),lower(tag)))
      fprintf('%s: Skipping unknown tag: %s\n', mfilename, tag);
      continue
    end

    % Assign the field the value
    st.(tag) = args{iarg+1};
  end % for iarg
end % if nargin
