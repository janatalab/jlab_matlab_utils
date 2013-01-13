function [strs,cnt] = cellhist(vals)
% [strs,cnt] = cellhist(vals);
%
% Counts the number of instances of each string in the cell array vals
%

% 8/27/03 PJ

sorted = sort(vals);
unique_vals = unique(sorted);
nunique = length(unique_vals);

cnt = zeros(1,nunique);

for ival = 1:nunique
  cnt(ival) = sum(ismember(sorted,unique_vals{ival}));
end

strs = unique_vals;

return
