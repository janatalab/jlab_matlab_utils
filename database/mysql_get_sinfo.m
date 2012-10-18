function sinfo = mysql_get_sinfo(subid,params)

% Returns subject info for given subject id.
% 
%   sinfo = mysql_get_sinfo(subid,conn_id);
%
% Returns information about the subject from the database and puts it into a
% structure
%
% REQUIRES
%   subid - a subject ID for which to retrieve subject table info
%   params
%       .conn_id - active MySQL connection ID to use.
%       .enc_key - output of mysql_login, to decrypt subject table data
% 
% RETURNS
%   sinfo
%       .datenum -- date of birth
%       .gender
%       .exp_ids -- experiments the subject has participated in
%       .exp_names
%       .sess_ids -- session IDs

% 03/05/06 Petr Janata - determines participation dates for the subject
% 04/15/06 PJ - optimized search for subject ID in response tables
% 06/21/06 PJ - searches session table instead of response tables.
% 10/05/09 Stefan Tomic - reads in subject name and DOB as encrypted data
% 10/30/09 PJ - minor fix to convert date string to datenum, following
%               subject table encryption modification
% 06/15/10 PJ - mysql_make_conn
% 02/17/11 FB - adapted for new mysql authentication scheme, using
% mysql_login before this function to get login information for the
% particular user that will be running this function. This will
% 03/29/11 FB - if DOB is encrypted, datenum(sinfo.datenum) will fail. This
% error is now being caught, so that mysql_get_sinfo can continue
% 10/18/12 PJ - ignores experiments with no response table

sinfo = [];

% Do some parameter checking
msg = nargchk(1,2,nargin);
if ~isempty(msg)
  disp(msg)
  return
end

% Check for valid connection to database
if isempty(params) || ~isfield(params,'conn_id') || ...
        isempty(params.conn_id) || mysql(params.conn_id,'status')
  error('%s: Do not have a valid connection ID', mfilename);
end

% Initialize the output structure
sinfo = ...
    struct('datenum',[],'gender','','name_last','','name_first','','exp_ids',[],'sess_ids',[],'exp_date',[],'exp_names',{''});

if isfield(params,'enc_key') && ~isempty(params.enc_key)
  mysql_str = sprintf(['select date_entered, aes_decrypt(`dob`,''%s''),'...
      ' gender, aes_decrypt(`name_last`,''%s''), '...
      'aes_decrypt(`name_first`,''%s'')  from subject '...
      'where subject_id="%s";'],...
      params.enc_key,params.enc_key,params.enc_key,subid);
else
  mysql_str = sprintf(['select date_entered, dob, gender, name_last, '...
      'name_first from subject where subject_id="%s";'],subid);
end
[date_entered...
    sinfo.datenum, sinfo.gender, sinfo.name_last, sinfo.name_first] = ...
    mysql(params.conn_id, mysql_str);

% Convert date string to datenum
try sinfo.datenum = datenum(sinfo.datenum);
catch sinfo.datenum = 'encrypted'; end

%
% Determine which experiments this subject has been in.  Do this on the basis
% of the session table. Prior to 11/21/05, the subject ID was not being written
% into the session table, so for subjects who entered the database, we have to
% do things the slow way just to be safe.
%

crit_datenum = 732637; % datenum('Nov-21-2005');

if date_entered >= crit_datenum
  % Get a list of the experiment and session IDs as well as dates
  mysql_str = sprintf('SELECT experiment_id, session_id, date_time  FROM session WHERE subject_id="%s";', subid);
  [sinfo.exp_ids, sinfo.sess_ids, sinfo.exp_date] = mysql(params.conn_id, mysql_str);
  
  % Get the experiment names
  expids = sprintf('%d,', sinfo.exp_ids(:));
  expids(end) = [];
  mysql_str = sprintf(['SELECT experiment_title FROM experiment WHERE' ...
	' experiment_id IN (%s);'], expids);
  [sinfo.exp_names] = mysql(params.conn_id, mysql_str);
  
else

  % First, get the list of response tables
  mysql_str = 'SELECT experiment_id, experiment_title, response_table FROM experiment;';
  [expids, exptitles, resp_tbls] = mysql(params.conn_id, mysql_str);

  % Now, loop through the tables, checking whether this particular subject has
  % entries in this table
  ntbl = length(resp_tbls);
  tbl_mask = zeros(1,ntbl);
  for itbl = 1:ntbl
		if isempty(resp_tbls{itbl})
			continue
		end
    mysql_str = sprintf('SELECT subject_id FROM %s WHERE subject_id="%s";', resp_tbls{itbl}, subid);
    sublist = mysql(params.conn_id, mysql_str);
    if ~isempty(sublist)
      tbl_mask(itbl) = 1;
    end
  end % for itbl = 1:ntbl

  exp_idxs = find(tbl_mask);
  sinfo.exp_ids = expids(exp_idxs);
  sinfo.exp_names = exptitles(exp_idxs);

  % Get participation dates
  for itbl = 1:length(exp_idxs)
    mysql_str = sprintf('SELECT date_time, session_id FROM %s where subject_id="%s"', resp_tbls{exp_idxs(itbl)}, subid);
    [dateinfo, session_ids] = mysql(params.conn_id, mysql_str);
    [sessids,sess_offsets] = unique(session_ids);
    for isess = 1:length(sessids)
      sinfo.exp_date{itbl}(isess) = dateinfo(sess_offsets(isess));
    end
  end % for itbl = 1:length(exp_idxs
  
end % if date_entered >= crit_datenum
