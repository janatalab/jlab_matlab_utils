function [user,passwd,enc_key] = ensemble_get_credentials(auth_path,cred_type)
% [user,passwd,enc_key] = ensemble_get_credentials(auth_path);
%
% Looks for a file called <cred_type>_database_credentials.txt in the directory
% specified in authpath, where cred_type is either researcher or subject

% 06/15/10 PJ - fixed original version

user = '';
passwd = '';
enc_key = '';

if nargin < 1
  auth_path = './';
end

if ~exist('cred_type','var')
  cred_type = 'subject';
elseif ~ismember(cred_type, {'subject','researcher'})
  fprintf('%s: Unknown login type: %s\n', mfilename, cred_type);
  return
end

% Make sure the directory exists
if ~exist(auth_path,'dir')
  fprintf('%s: Directory does not exist: %s\n', mfilename, auth_path);
  return
end

% Make sure the credentials file exists and is readable
fname = fullfile(auth_path,sprintf('%s_database_credentials.txt', cred_type));
if ~exist(fname,'file')
  fprintf('%s: Credentials file does not exist or is not readable.\n', mfilename);
  return
end

fid = fopen(fname,'rt');
if fid==-1
  warning(sprintf('Failed to open credentials file: %s\n', fname));
  return
end

while ~feof(fid)
  str = fgetl(fid);
  [tag,val] = strtok(str,':');
  
  switch tag
    case 'user'
      user = strtrim(val(2:end));
    
    case 'passwd'
      passwd = strtrim(val(2:end));
      
    case 'enc_key_file'
      fstub = strtrim(val(2:end));
      if isempty(fstub)
        continue
      end
      enc_fname = fullfile(auth_path,fstub);
      encfid = fopen(enc_fname,'rt');
      if encfid < 3
        enc_key = '';
        fprintf('Could not open encryption key file\n');
        continue
      else
        enc_key = fgetl(encfid);
        fclose(encfid);
      end
      
    end
  
end

fclose(fid);

