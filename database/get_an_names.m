function names = get_an_names(an,verbose)
% Returns list of analysis names in an an{na} array
%
% names = get_an_names(an,verbose);

% January 5, 2011 PJ
% 05Mar2014 PJ - added verbosity option

if nargin < 1
  fprintf('Need to pass in an analysis stack\n')
  return
end

if nargin < 2
  verbose = true;
end

na = length(an);
names = cell(na,1);
for ia = 1:na
	names{ia} = an{ia}.name;
  if verbose || nargout < 1
    fprintf('%d: %s\n', ia, names{ia});
  end
end

return
