function [str] = prob2str(prob,crit,tok)
% [str] = prob2str(prob,crit,tok)
%
% Converts probability value in prob to a sequence of tokens.
%
% E.g. 0.001 -> ***
%

% 03/05/05/ PJ

if nargin < 3
  tok = '*';
end

if prob < 1E-5
  str = repmat(tok,1,5);
elseif prob <= crit
  str = repmat(tok,1,abs(fix(log10(prob))));
else
  str = 'n.s.';
end