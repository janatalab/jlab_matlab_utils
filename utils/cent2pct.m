function pct = cent2pct(cent)
% cent2pct.m
%
% Converts cents to percent
%
% Plots percent deviation as a function of cents deviation if a vector of cents
% values is given
%
% 1 semitone = 100 cents

% 05/20/09 Petr Janata

pct = (2.^(cent/1200)-1)*100;

if length(cent) > 1
  plot(cent,pct)
  xlabel('Cents')
  ylabel('Percent')
end
