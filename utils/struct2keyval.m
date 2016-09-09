function outcell = struct2keyval(instruct)
% Converts a structure into a cell array that can be used as a series of
% key/value pairs that can be passed into functions, i.e. as varargin
%
% outcell = struct2keyval(instruct);
%
% outcell is a 1 X N cell array where N is twice the number of fields in
% the structure instruct.

% 09Sep2016 Petr Janata

outcell = [fieldnames(instruct)'; struct2cell(instruct)'];
outcell = outcell(:)';

end
