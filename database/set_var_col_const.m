function col = set_var_col_const(var_list)
% Assigns an index to a field that is named after given variable
%
% col = set_var_col_const(var_list);
%
% This is a generic mechanism for assigning an index to a field that is named
% after a variable in the variable list. This allows for scripts to index into
% data using names and saves the script author from having to do strmatch()
% lookups each time.
%
% See also: set_form_col_const.m for a fixed set of column IDs
%
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved
%
% Author:
% 01/26/07 Petr Janata

for ivar = 1:length(var_list)
  col.(var_list{ivar}) = ivar;
end
