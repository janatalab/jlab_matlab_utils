function names = get_an_names(an)
% Returns list of analysis names in an an{na} array
%
% names = get_an_names(an);

% January 5, 2011 PJ

if nargin < 1
  fprintf('Need to pass in an analysis stack\n')
  return
end

na = length(an);
names = cell(na,1);
for ia = 1:na
	names{ia} = an{ia}.name;
	fprintf('%d: %s\n', ia, names{ia});
end

return
