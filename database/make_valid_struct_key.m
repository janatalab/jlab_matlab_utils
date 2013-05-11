function [str] = make_valid_struct_key(str)

% [str] = make_valid_struct_key(str)
% 
% given str, will replace all non alpha-numeric characters with underscores
% and will append 's' to the beginning of the string
% 
% this will construct valid struct field names from any string
% 
% NOTE: if you have two strings that are only differentiated by their
% non alpha-numeric characters, this will not maintain their unique
% characters
% 

% fb 2007/10/29

% if (~isstr(str))
%     error(sprintf('the variable that you passed is not a string'));
% end

str = regexprep(str,'\"','');
str = regexprep(str,'[^a-zA-Z0-9]','_');
str = strcat('s',str);
