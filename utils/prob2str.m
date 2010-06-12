function [str] = prob2str(prob,crit,tok,places)
% [str] = prob2str(prob,crit,tok,places)
%
% Converts probability value in prob to a sequence of tokens.
%
% crit - criterion probability above which to return 'n.s.'
% tok - token to use for string. Default: *
% places - number of places to go out to. Default: 1E-4 = ****
%
% E.g. 0.001 -> ***
%

% 03/05/05/ PJ

if nargin < 4
  places = 4;
end

if nargin < 3
  tok = '*';
end

if prob < 10^-places
  str = repmat(tok,1,places);
elseif prob <= crit
  str = repmat(tok,1,abs(fix(log10(prob))));
else
  str = 'n.s.';
end