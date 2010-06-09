% ensemble_globals_template.m
%
% Janata Lab, Center for Mind and Brain, UC Davis
%
% Original Version by Petr Janata
% Edited June 8, 2010 Stefan Tomic - created template from original file. Added 'county' to list of encrypted fields
%
% Copy this file to a file called "ensemble_globals.m"
% Then edit the path variables to point to locations relevant for your Ensemble environment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill in the following paths and variables

%the default hostname and database for Ensemble database connections
%if host and/or database are omitted when calling mysql_make_conn,
%the values will automatically be populated with these defaults.
DEFAULT_HOST = '';
DEFAULT_DATABASE = '';

%path to the root directory of all ensemble stimuli (e.g. /var/www/html/ensemble/stimuli)
stimulus_root = ''; 

%Path to the ipem analysis directory (which mimics the structure of stimulus_root. 
%You do not need this if you do not use our scripts that call the IPEM Toolbox
stimulus_ipem_analysis_root = ''; 

%Path to the Ensemble encryption key file (for encrypting sensitive subject data)
paths.data_encryption_key = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paths.stimulus_root = stimulus_root;
paths.stimulus_ipem_analysis_root = stimulus_ipem_analysis_root;

encrypted_fields.subject = {'name_last','name_first','name_middle','name_suffix','email','phone1','phone2',...
		    'address1','address2','address3','city','county','state','postal_code','dob'};
