function params = ensemble_check_param(checkParam,params,default)
% Checks for the existance of checkParam within params.
%
% ensemble_check_param(checkParam,params,default)
%
% sees if checkParam, passed in as a string, is set in the params
% structure. Each param structure 'level' will be traversed and
% checked for existence.
%
% For example, if checkParam = 'params.test1.test2.test3' and
% params.test1 exists but not the other fields, this function will
% traverse each level (e.g. params.test1.test2 = [], then
% params.test1.test2.test3 = []. Only fields that don't exist yet are
% initialized and set to null.
%
% If a default value is given and the field wasn't already assigned
% a value, then the default value will be assigned to the
% last level (e.g. params.test1.test2.test3 = 4 if the default was
% 4)
%
% if one of the lower params levels exists and is equal to a value, then
% the function will stop processing and return an error, instead of
% overwriting that value with a field (e.g. if params.test1 = 2,
% then the function will not attempt to overwrite it with
% params.test1.test2 = [])
%
% This function was primarily motivated by the desire to avoid lots
% of if-then-else statements when checking for parameter
% default values in ensemble functions.
%
% 30 March 2007 - Stefan Tomic



parm = regexp(checkParam,'\.(\w*[^\.])','tokens');
nLevels = length(parm);

%the string that will be used to check field existence
compare = 'params';

%iLevel refers to each param struct level
for iLevel = 1:nLevels

  %see if the field is set for the next level
  eval(sprintf('paramExists = isfield(%s,parm{%d}{1});',compare,iLevel));

  %also see if the current level contains a value, in which case we
  %can't overwrite it with the field for the next level
  eval(sprintf('valueSet = ~isempty(%s) & ~isstruct(%s);', ...
	       compare,compare));
    
  %throw an error if the current level contains a value
  if(valueSet)
    error(sprintf('%s contains a value so parameter checking cannot proceed',compare));
  end

  
  %if the field for the next level does not exist, then create it
  %and set it to empty
  if(~paramExists)
    eval(sprintf('%s.(parm{%d}{1}) = [];',compare,iLevel));
  end
  
  %append the field for the next level to the compare string
  compare = [compare sprintf('.%s',parm{iLevel}{1})];  
  
end

%now that we are at the last level, check if what we think should
%be the last level actually contains a struct. If so, throw a
%warning (since it's not empty, it will not be overwritten by the
%following code)
eval(sprintf('setToStruct = isstruct(%s);',compare));
if(setToStruct)
  warning(sprintf('%s contains a structure. Parameter checking cannot proceed'));
end

%finally see if a value is already set at the last level
eval(sprintf('valueSet = ~isempty(%s) & ~setToStruct;',compare));

%if a value is not set at the last level and a default value has
%been submitted, then set it to the default value
if(exist('default','var') & ~valueSet)
  eval(sprintf('%s = default;',compare));
end
