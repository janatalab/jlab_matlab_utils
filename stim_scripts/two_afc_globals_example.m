function params = two_afc_globals
% params = two_afc_globals;
%
% Copyright (c) 2006-13 The Regents of the University of California, Davis campus. All Rights Reserved.
%
% Author: Petr Janata
% Edits by JR

params.version = 'two_afc_globals';
params.debug = 3;

params.mysql.host = 'localhost';      % replace with database host name
% Database name should match an entry in the mysql_login.m file.
params.mysql.database = 'ensemble';   % replace with database name
params.mysql.conn_id = 3;

params.mysql.login_type = 'subject';

params.mysql = mysql_login(params.mysql);  % login to the database

% list of relevant stimulus attributes in the database table 'trial_x_attribute'
% used to look up the trials to include in this experiment
params.attrib_name = '2afc_example';

% boolean: determine whether to show all trials or just a subset
% Current two_afc_select.m only supports value of true for this parameter
params.include_all_trials = 1;

% variables for analysis and logging
params.ensemble.expname = '2afc_example';      % experiment name as defined in the QEI
params.resptbl_name = 'response_2afc_example'; % response table name defined in QEI
