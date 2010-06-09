% ensemble_researcher_login_template.m
%
% Janata Lab, Center for Mind and Brain, UC Davis
%
% Original Version by Petr Janata
% Edited June 8, 2010 Stefan Tomic - created template from original file. Added 'county' to list of encrypted fields
%
% Copy this file to a file called "mysql_researcher_login.m"
% and save it in a location that is only readable by users that run Ensemble analyses
% This file should not be readable by others.
% Change <hostname> to the hostname where your MySQL database resides (e.g. localhost).
% Change <database> to the name of your ensemble database
% Change <username> to your primary Ensemble administration account (e.g. experimenter)
% Change <password> to the password that corresponds to this account.
%
% The switch statements provide for conditional branching based on hostname and
% database. If you are supporting database connections to
% multiple hostnames and databases, different user and passwd values can be
% populated for them.


switch host
  case 	{'<hostname>'}
   switch db
    case {'<database>'}
     user = '<username>';
     passwd = '<password>';
   end
end