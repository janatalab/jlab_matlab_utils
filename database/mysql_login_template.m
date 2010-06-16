% function params=mysql_login(params)
%
% ensemble_researcher_login_template.m
%
% Janata Lab, Center for Mind and Brain, UC Davis
%
% Original Version by Petr Janata
% Edited June 8, 2010 Stefan Tomic - created template from original file. Added 'county' to list of encrypted fields
%
% Copy this file to a file called "mysql_login.m"
% and save it in a location that is only readable by users that run Ensemble analyses
% This file should not be readable by others.
%
% params is a structure containing the following fields:
% .host to the hostname where your MySQL database resides (e.g. localhost).
% .database to the name of your ensemble database
%
% The switch statements provide for conditional branching based on hostname and
% database. If you are supporting database connections to
% multiple hostnames and databases, different user and passwd values can be
% populated for them.

%
% Add code to check integrity of .host and .database fields and populate
% with defaults if necessary.
%

% End of code that checks for host and databse names

% Check to see if we are going to authenticate with subject or researcher
% privileges. Default to more restrictive subject privileges
if ~isfield(params,'login_type') || isempty(params.login_type)
  params.login_type = 'subject'; 
end

% Determine where files with authentication information live
try auth_root = params.auth_root; catch params.auth_root = '/var/private/'; end
switch params.host
  case 	{'127.0.0.1','localhost'}
   switch db
     case {'ensemble_main'}
       auth_path = fullfile(params.auth_root,'main');

     case {'ensemble_dev','ensemble_test'}
       auth_path = fullfile(params.auth_root,'dev');
   
   end % switch db

end % switch host

% Try to read the authentication information
[user, passwd, enc_key] = ensemble_get_credentials(auth_path,params.login_type);

% Copy to output structure
params.user = user;
params.passwd = passwd;
params.enc_key = enc_key;

return
