function fh = parse_fh(target)

% returns a function handle from a target (string, function library, or function itself)
% 
%   fh = parse_fh(target)
% 
% if given a string, uses str2func to retrieve function handle
% if given a function handle, returns that function handle
% 
% function library: a function that, given an argument, will return the
% function handle for a sub-function of that function.
% 
% if given a cell array of 2 strings, parse_fh  will attempt to turn the
% first string into a function handle (str2func) and then will pass the
% second string as an argument to the first string/function, and will
% expect a function handle to return as the first product of this function
% this is the convention that is being used in stimulus selection function
% libraries (see private/matlab/projects/stim_markup/stim_markup_library.m
% for an example)
% 
% FB 2009.04.20

if ischar(target)
  % target is a string, use str2func
  fh = str2func(target);
elseif iscell(target)
  % target is a cell ... could this be intended for a function library?
  if length(target) ~= 2 || ...
          ~ischar(target{1}) || ...
          ~ischar(target{2})
    error(['must provide library name and sub-function for '...
        'function library reference\n']);
  end
  lfh = str2func(target{1});
  fh = lfh(target{2});
elseif isa(target,'function_handle')
  fh = target;
else
  error(['must provide either a function string, cell array of strings '...
      'for a function library call, or a valid function\n']);
end
