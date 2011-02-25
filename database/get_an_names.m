function names = get_an_names(an)
% Returns list of analysis names in an an{na} array
%
% names = get_an_names(an);

% January 5, 2011 PJ

na = length(an);
names = cell(na,1);
for ia = 1:na
	names{ia} = an{ia}.name;
end

return