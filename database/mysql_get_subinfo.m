function subinfo = mysql_get_subinfo(varargin)
% Returns subject info for given subjects
% subinfo = mysql_get_subinfo(varargin)
%
% Returns information associated with one or more subjects provided in an array
% of subjects. Arguments are passed in as tag/value pairs.
%
% In order to get information from encrypted fields in the subject table it
% is necessary to pass in encryption key information. The encryption key is
% obtained using mysql_researcher_login.
%
% Supported input arguments (tags):
% 'subject_id' - vector of subject IDs
% 'conn_id' - mysql connection ID to utilize
% 'mysql' - a structure containing host,database,user,passwd,enc_key fields
% 'extract_vars' - list of fields to extract

% 01/26/07 Petr Janata
% 06/01/07 Stefan Tomic - fixed handling of temporary conn_id
% 06/15/10 PJ - mysql_make_conn handling

ensemble_globals;
  
% Initialize some variables
subinfo = ensemble_init_data_struct;

extract_vars = {};

% Parse the input arguments
narg = length(varargin);
for iarg = 1:2:narg
  switch varargin{iarg}
    case {'subject_id','subject_ids'}
      subids = varargin{iarg+1};
    case {'mysql','ensemble'}
      params = varargin{iarg+1};
		case {'extract_vars'}
			extract_vars = varargin{iarg+1};
    otherwise
      fprintf('mysql_get_subinfo: Unknown input argument: %s\n', varargin{iarg});
  end
end

if isfield(params,'conn_id')
	conn_id = params.conn_id;
else 
	conn_id = [];
end

if isfield(params,'enc_key')
	enc_key = params.enc_key;
else 
	enc_key = '';
end

% Try to get encryption key if this hasn't been obtained yet
tmp_conn_id = false;
if isempty(conn_id) || isempty(enc_key)
	if ~exist('params','var') || ~isfield(params,'conn_id') || ...
			~isfield(params,'user') || ~isfield(params,'passwd') || ...
			~isfield(params,'host') || ~isfield(params,'database')
		error('%s: Insufficient information for establishing database connection', mfilename)
	end
	
	params.login_type = 'researcher';
	params = mysql_login(params);
	
	% Check for connection to database
	
	if mysql(params.conn_id,'status')
		tmp_conn_id = true;
  else
    tmp_conn_id = false;
	end
	
	conn_id = mysql_make_conn(params);
end

subinfo.type = 'subject_info';

% Call mysql_extract_data
enc_fields = encrypted_fields.subject;

[subinfo.data, subinfo.vars] = mysql_extract_data('table','subject',...
  'subject_id',subids,...
  'encrypted_fields',enc_fields,...
  'enc_key',params.enc_key, ...
	'extract_vars', extract_vars, ...
  'conn_id',conn_id);

% Close the database connection if it was temporary
if tmp_conn_id
  mysql(conn_id,'close');
end

return
