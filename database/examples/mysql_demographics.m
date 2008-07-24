function [dg, subject_vars] = mysql_demographics(params, conn_id)
% [dg, subject_vars] = mysql_demographics(params, conn_id);
%
% Gets subject demographic information based on the information request
% specified in the params structure
%
% params
%    .exp_list -- cell string of experiments to pull data for
%    .irb_id -- protocol numbers
%    .extract_vars -- list of variables to extract
%    .exclude_subs -- list of subject IDs to exclude from the summaries,
%                     e.g. test subject IDs, etc.
%    .start_date -- start of reporting period (must be in datenum format)
%    .stop_date -- end of reporting period (must be in datenum format)
%
% conn_id -- mysql connection ID
%

% 01/25/05 Petr Janata
% 06/26/05 PJ -- cleaned up conn_id handling

CONN_ID = 1;

dg = {};

% Data format table identifiers for enums containing the different category labels
GENDER_DFID = 77;
ETHNICITY_DFID = 78;
RACE_DFID = 79;

GENDER_QID = 317;
ETHNICITY_QID = 318;
RACE_QID = 319;
DEMO_VAR_LIST = {'gender','ethnicity','race'};

race_mapping = {...
      'Caucasian', {'Cauc','Caucasian'}; ...
      'AfrAmer', {'Black or African American', 'AfrAm'}; ...
      'Asian', {'Asian'}; ...
      'NatAmer', {'American Indian','NatAm'}; ...
      'Hawaiian', {'Native Hawaiian','NatHaw'}; ...
      'Multiple', {'More','More than One Race'}; ...
      'Other', {'Other'}; ...
      'Unknown', {'Unknown','Unknown or No Answer'}};


% Check the input arguments
if isempty(params.exp_list)
  exp_list = {};
  INFO_BY_EXP = 0;
else
  exp_list = params.exp_list;
  INFO_BY_EXP = 1;
end

if isempty(params.irb_id)
  irb_list = {};
  INFO_BY_IRB = 0;
else
  irb_list = params.irb_id;
  INFO_BY_IRB = 1;
end

if ~isfield(params,'exclude_subs') || isempty(params.exclude_subs)
  exclude_subs = {};
else
  exclude_subs = params.exclude_subs;
end

if ~isfield(params,'start_date') || isempty(params.start_date)
  start_date = datenum('01-Jul-2004');
else
  start_date = params.start_date;
end

if ~isfield(params,'stop_date') || isempty(params.stop_date)
  stop_date = datenum(now);
else
  stop_date = params.stop_date;
end


% List of subject variables to extract
subject_vars = { ...
      'subject_id'; ...
      'gender'; ...
      'ethnic_category'; ...
      'racial_category'; ...
      'dob'; ...
      };

% Do some additional variable checking
if isstr(exp_list)
  exp_list = {exp_list};
end

if isstr(irb_list)
  irb_list = {irb_list};
end

% Connect to default host and database if no connection ID is specified
try conn_id(1);
catch 
  tmp_conn_id = 1;
  mysql_make_conn([],[],CONN_ID);
  conn_id = CONN_ID;
end

% Initialize some variables
exp_id = [];
exp_names = {};
resp_tbl_list = {};

% Get a list of the response tables we need to deal with
if INFO_BY_EXP
  sql_str = sprintf(['SELECT experiment_id, experiment_title, response_table FROM experiment' ...
	' WHERE experiment_title IN ("%s");'], cell2str(exp_list,'","'));
  [tmp_id,tmp_names,tmp_tbl_list] = mysql(conn_id,sql_str);

  fprintf('Found %d/%d experiments matching the experiment list:\n%s\n', length(tmp_names), length(exp_list), cell2str(tmp_names,'\n'));
  
  exp_id = [exp_id, tmp_id];
  exp_names = [exp_names, tmp_names];
  resp_tbl_list = [resp_tbl_list, tmp_tbl_list];
end % if INFO_BY_EXP

if INFO_BY_IRB
  sql_str = sprintf(['SELECT experiment_id, experiment_title, response_table FROM experiment' ...
	' WHERE irb_id IN (%s);'], cell2str(irb_list,','));
  [tmp_id,tmp_names,tmp_tbl_list] = mysql(conn_id,sql_str);
  nexp = length(tmp_names);

  fprintf('Found %d experiments matching the IRB ID list (%s):\n', length(tmp_names), cell2str(irb_list,','));
  for iexp = 1:nexp
    fprintf('%s\t%d\n',tmp_names{iexp}, tmp_id(iexp));
  end
  
  exp_id = [exp_id, tmp_id];
  exp_names = [exp_names, tmp_names];
  resp_tbl_list = [resp_tbl_list, tmp_tbl_list];
end % if INFO_BY_IRB

% In case we were selecting by both experiment title and IRB protocol ID, prune
% the list for only the unique experiments
[unique_ids, unique_idxs] = unique(exp_id);
exp_id = exp_id(unique_idxs);
exp_names = exp_names(unique_idxs);
resp_tbl_list = resp_tbl_list(unique_idxs);

% Now generate a list of subject IDs
exp_subs = {};
sub_exp_ids = [];
sess_ids = [];
for iexp = 1:length(exp_names)
  % Get all the subject IDs and session IDs in this response table
  [tmp_subid, tmp_sessid] = mysql(conn_id,sprintf('SELECT subject_id, session_id FROM %s WHERE experiment_id=%d;', resp_tbl_list{iexp}, exp_id(iexp)));
  
  % Match on unique session_id 
  [tmp_sessid, tmp_idxs] = unique(tmp_sessid);
  sess_ids = [sess_ids; tmp_sessid];
  exp_subs = [exp_subs; tmp_subid(tmp_idxs)];
  sub_exp_ids = [sub_exp_ids; ones(size(tmp_sessid))*exp_id(iexp)];
end

rm_idxs = [];
% Strip out any temporary subject IDs
rm_idxs = [rm_idxs; strmatch('tmp_',exp_subs)];

% Remove any unwanted subject IDs
if ~isempty(exclude_subs)
  [rm_mask] = ismember(exp_subs,exclude_subs);
  rm_idxs = [rm_idxs; find(rm_mask)];
end
exp_subs(rm_idxs) = [];
sub_exp_ids(rm_idxs) = [];
sess_ids(rm_idxs) = [];

% Generate a comma-separated list of subject IDs
sub_str = sprintf('''%s'',', exp_subs{:});
sub_str(end) = [];

% Extract information from the subject table. Note, following the revision of
% the handling of demographics forms, demographic information is stored in the
% response tables, so we really have to look in two place to get all the info.
sql_str = sprintf('SELECT %s FROM subject WHERE subject_id IN (%s)', ...
    cell2str(subject_vars,','), sub_str);

% Get the data from the database
[dg{1:length(subject_vars)}] = mysql(conn_id,sql_str);

% Figure out which of the subjects we need to look into the response tables for
% the demographic info
tmp_sub_ids = dg{:,strcmp(subject_vars,'subject_id')};
tmp_gender = dg{:,strcmp(subject_vars,'gender')};
tmp_ethnicity = dg{:,strcmp(subject_vars,'ethnic_category')};
tmp_race = dg{:,strcmp(subject_vars,'racial_category')};

gender_empty_mask = cellfun('isempty',tmp_gender);
ethnicity_empty_mask = cellfun('isempty', tmp_ethnicity);
race_empty_mask = cellfun('isempty', tmp_race);

% create the missing data mask
sub_tbl_empty_mask = gender_empty_mask & ethnicity_empty_mask & ...
    race_empty_mask;

missing_subs_mask = ismember(exp_subs, tmp_sub_ids(sub_tbl_empty_mask));
missing_subs = exp_subs(missing_subs_mask);
missing_subs_str = sprintf('''%s'',', missing_subs{:});
missing_subs_str(end) = [];

% Remove the entries for these subjects from the current data
for idg = 1:length(dg)
  dg{idg}(sub_tbl_empty_mask,:) = [];
end

sub_ids = exp_subs(~missing_subs_mask);
sess_ids = sess_ids(~missing_subs_mask);

% get list of experiments whose response tables we need to check
new_exp_id_list = unique(sub_exp_ids(missing_subs_mask));
sub_exp_ids = sub_exp_ids(~missing_subs_mask);

% Now extract variables from the data structure for existing subjects
tmp_sub_ids = dg{:,strcmp(subject_vars,'subject_id')};

nsub_ids = length(sub_ids);
gender = cell(nsub_ids,1);
ethnicity = cell(nsub_ids,1);
race = cell(nsub_ids,1);
age = zeros(nsub_ids,1);
unique_subids = unique(sub_ids);
nunique = length(unique_subids);
for isub = 1:nunique
  curr_subid = unique_subids(isub);
  curr_sub_ids_idxs = strmatch(curr_subid,sub_ids,'exact');
  curr_dg_idx = strmatch(curr_subid,tmp_sub_ids,'exact');
  gender(curr_sub_ids_idxs) = dg{:,strcmp(subject_vars,'gender')}(curr_dg_idx);
  ethnicity(curr_sub_ids_idxs) = dg{:,strcmp(subject_vars,'ethnic_category')}(curr_dg_idx);
  race(curr_sub_ids_idxs) = dg{:,strcmp(subject_vars,'racial_category')}(curr_dg_idx);
  age(curr_sub_ids_idxs) = dg{:,strcmp(subject_vars,'dob')}(curr_dg_idx);
end

% Retrieve the category labels
[cat_labels] = mysql_resolve_enum([GENDER_DFID, ETHNICITY_DFID, RACE_DFID], conn_id);

% Now loop over the experiments that we need to consult in order to pull in the
% demographic information from the response tables
new_sub_ids = {};
new_sess_ids = [];
for iexp = 1:length(new_exp_id_list)
  curr_exp_id = new_exp_id_list(iexp);
  curr_resp_tbl = resp_tbl_list{find(exp_id == curr_exp_id)};
  
  sql_str = sprintf(['SELECT subject_id, response_enum, question_id, session_id FROM %s ' ...
	'WHERE experiment_id=%d AND question_id IN (%d,%d,%d) ' ...
	'AND subject_id IN (%s);'], curr_resp_tbl, curr_exp_id, GENDER_QID, ETHNICITY_QID, RACE_QID, missing_subs_str);
  [tmp_subid, tmp_resp, tmp_qid, tmp_sessid] = mysql(conn_id,sql_str);
  
  new_sub_ids = [new_sub_ids; tmp_subid(1:3:end)];
  new_sess_ids = [new_sess_ids; tmp_sessid(1:3:end)];
  
  for ivar = 1:length(DEMO_VAR_LIST)
    curr_var = DEMO_VAR_LIST{ivar};
    tmp_mask = [];
    switch curr_var
      case 'gender'
	tmp_mask = (tmp_qid == GENDER_QID);
      case 'ethnicity'
	tmp_mask = (tmp_qid == ETHNICITY_QID);
      case 'race'
	tmp_mask = (tmp_qid == RACE_QID);
    end
    
    % Populate the gender, ethnicity, and race variables
    enum_idxs = log2(tmp_resp(tmp_mask))+1;
    tmp_vals = cat_labels{ivar}(enum_idxs);
    tmp_vals = strrep(tmp_vals,'"','');
    eval([curr_var '=[' curr_var '; tmp_vals''];']);
    
  end % for ivar=
end % for iexp = 1:length(new_exp_id_list)

sub_ids = [sub_ids; new_sub_ids];
sess_ids = [sess_ids; new_sess_ids];

% Get the date of participation from the session table
sessid_str = sprintf('%d,', sess_ids);
sessid_str(end) = [];
sql_str = sprintf('SELECT date_time FROM session WHERE session_id IN (%s);', sessid_str);
session_date = mysql(conn_id,sql_str);

date_mask = (session_date >= start_date) & (session_date <= stop_date);
valid_mask = date_mask;

%
% Print some summary info
%
fprintf('\nTotal number of subjects: %d\n', sum(valid_mask));
fprintf('Start date: %s\n', datestr(start_date));
fprintf('Stop date: %s\n\n', datestr(stop_date));

male_mask = ismember(gender,{'M','Male'}) & valid_mask;
female_mask = ismember(gender,{'F','Female'}) & valid_mask;
unknown_gender_mask = ~male_mask & ~female_mask & valid_mask;

% Hispanic breakdown
hispanic_mask = ismember(ethnicity,{'HL','Hispanic or Latino'}) & valid_mask;
nonhispanic_mask = ismember(ethnicity,{'Not Hispanic or Latino','not_HL'}) & valid_mask;
hispanic_unkwn_mask = ~hispanic_mask & ~nonhispanic_mask & valid_mask;

fprintf('\t%15s\t%15s\t%15s\n','Hispanic','Non-Hispanic','Unknown');
fprintf('Male:\t%15d\t%15d\t%15d\n', ...
    sum(male_mask & hispanic_mask), ...
    sum(male_mask & nonhispanic_mask), ...
    sum(male_mask & hispanic_unkwn_mask));
fprintf('Female:\t%15d\t%15d\t%15d\n', ...
    sum(female_mask & hispanic_mask), ...
    sum(female_mask & nonhispanic_mask), ...
    sum(female_mask & hispanic_unkwn_mask));
fprintf('Unknown:%15d\t%15d\t%15d\n', ...
    sum(unknown_gender_mask & hispanic_mask), ...
    sum(unknown_gender_mask & nonhispanic_mask), ...
    sum(unknown_gender_mask & hispanic_unkwn_mask));
fprintf('Total:\t%15d\t%15d\t%15d\n', sum(hispanic_mask), sum(nonhispanic_mask), sum(hispanic_unkwn_mask));

% Do overall racial breakdown
race_cat = race_mapping(:,1);
header = sprintf('\n\n\t');
for ir = 1:length(race_cat)
  header = sprintf('%s\t%8s', header,race_cat{ir});
  race_mask = ismember(race,race_mapping{ir,2});
  males_by_race(ir) = sum(race_mask&male_mask & valid_mask);
  females_by_race(ir) = sum(race_mask&female_mask & valid_mask);
  unknown_by_race(ir) = sum(race_mask&unknown_gender_mask & valid_mask);
  
  hispanic_males_by_race(ir) = sum(race_mask & male_mask & valid_mask & hispanic_mask);
  hispanic_females_by_race(ir) = sum(race_mask&female_mask & valid_mask & hispanic_mask);
  hispanic_unknown_by_race(ir) = sum(race_mask&unknown_gender_mask & valid_mask & hispanic_mask);
end
header = sprintf('%s\t%8s', header, 'Total');
males_by_race(end+1) = sum(males_by_race);
females_by_race(end+1) = sum(females_by_race);
unknown_by_race(end+1) = sum(unknown_by_race);
hispanic_males_by_race(end+1) = sum(hispanic_males_by_race);
hispanic_females_by_race(end+1) = sum(hispanic_females_by_race);
hispanic_unknown_by_race(end+1) = sum(hispanic_unknown_by_race);

fprintf('%s\n', header);
fprintf('Male:   %s\n', sprintf('\t%8d',males_by_race));
fprintf('Female: %s\n', sprintf('\t%8d',females_by_race));
fprintf('Unknown:%s\n', sprintf('\t%8d',unknown_by_race));
fprintf('Total:  %s\n', sprintf('\t%8d',males_by_race+females_by_race+unknown_by_race))

% Do hispanic racial breakdown
fprintf('\n\nRacial breakdown for hispanics')
fprintf('%s\n', header);
fprintf('Male:   %s\n', sprintf('\t%8d',hispanic_males_by_race));
fprintf('Female: %s\n', sprintf('\t%8d',hispanic_females_by_race));
fprintf('Unknown:%s\n', sprintf('\t%8d',hispanic_unknown_by_race));
fprintf('Total:  %s\n', sprintf('\t%8d',hispanic_males_by_race+hispanic_females_by_race+hispanic_unknown_by_race))

% Close the mysql connection
if exist('tmp_conn_id','var')
  mysql('close');
end
