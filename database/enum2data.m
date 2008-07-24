function data=enum2data(enum)
% Converts Ensemble enum values to numeric values.
%
% data=enum2data(enum);
%
% Any values that were left empty (0) are filled with NaN.

% 07/30/06 Petr Janata

% Convert the enum values into scale values. 
% Suppress log of zero errors in the process

warning off
data = log2(enum)+1;
warning on

data(isinf(data)) = NaN;

return
