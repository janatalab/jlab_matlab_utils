function subinfo = mysql_get_subinfo(varargin)
% Returns subject info for given subjects
% subinfo = mysql_get_subinfo(varargin)
%
% Returns information associated with one or more subjects provided in an array
% of subjects. Arguments are passed in as tag/value pairs.
%
% Supported input arguments (tags):
% 'subject_id' - vector of subject IDs
% 'conn_id' - mysql connection ID to utilize

% 01/26/07 Petr Janata
% 06/01/07 Stefan Tomic - fixed handling of temporary conn_id

ensemble_globals;
  
% Initialize some variables
subinfo = ensemble_init_data_struct;

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'subject_id','subject_ids'}
      subids = varargin{iarg+1};
    case 'conn_id'
      conn_id = varargin{iarg+1};
    case 'enc_key'
      enc_key = varargin{iarg+1};
    otherwise
      fprintf('mysql_get_subinfo: Unknown input argument: %s\n', varargin{iarg});
  end
end

% Check for connection to database
try conn_id(1);
catch   
  mysql_make_conn;
  tmp_conn_id = 1;
  conn_id = 0;
end

subinfo.type = 'subject_info';

% Call mysql_extract_data
enc_fields = encrypted_fields.subject;
if ~exist('enc_key')
  enc_key = ensemble_get_encryption_key;
end

[subinfo.data, subinfo.vars] = mysql_extract_data('table','subject','subject_id',subids,'encrypted_fields',enc_fields,'enc_key',enc_key,'conn_id',conn_id);

% Close the database connection if it was temporary
if (exist('tmp_conn_id','var'))
  mysql(conn_id,'close');
end
