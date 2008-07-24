function [turn_idxs, turn_vals] = find_turnarounds(cvals)
% [turn_idxs] = find_turnarounds(traj);
%
% Finds direction reversals in the input vector, and returns the indices of the
% reversals.
%

% 07/27/05 Petr Janata

cidxs = 1:length(cvals);

% Go ahead and get the turn arounds.  These occur whenever
% the direction of the bin values changes sign
cdir = diff(cvals);

% Remove repeated values
rep_idxs = find(~cdir);
cdir(rep_idxs) = [];
cidxs(rep_idxs+1) = [];

% Get the pure direction -- irrespective of jump size
dirvect = sign(cdir);

% Get indices that consistute turnarounds
turn_idxs = cidxs(find([0; diff(dirvect)])+1)-1;
turn_vals = cvals(turn_idxs);

return