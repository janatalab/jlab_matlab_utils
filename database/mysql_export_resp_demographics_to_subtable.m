function mysql_export_resp_demographics_to_subtable(varargin)
%
% Updates subject table from the demographics responses for a given experiment
% This is useful for older experiments which submitted demographics information
% to the response table instead of the subject table. The function populates
% the subject table from the responses.
%
% mysql_export_resp_demographics_to_subtable(varargin)
%
% INPUTS are tag/value pairs:
%
% 'exp_name' - the title of the experiment in the experiment table
% 'gender_qid' - the question ID of the gender question
% 'race_qid' - the question ID of the race question
% 'ethnicity_qid' - the question ID of the ethnicity question
% 'db' - the name of the database (if ommitted uses the default database)
% 'dryrun' - (optional) 0 or 1. If 1, performs a dry-run and just displays the sql
%                          commands to be executed rather than performing the
%                          updates. The default is 0.
%
%  Note that this function attempts to perform enum mapping between the
%  data_format table and the demographics fields in the subject table.
%  It is highly recommended that you perform a dryrun to ensure that the 
%  mapping is working (the mapping is hard-coded here, if there are other
%  maps they need to be introduced).
%
%  November 4, 2009 - Stefan Tomic, First Version
  
  
for iarg = 1:nargin
  
  switch(varargin{iarg})
    
   case 'exp_name'
    exp_name = varargin{iarg+1};
   case 'gender_qid'
    gID = varargin{iarg+1};
   case 'race_qid'
    rID = varargin{iarg+1};
   case 'ethnicity_qid'
    eID = varargin{iarg+1};
   case 'db'
    db = varargin{iarg+1};
   case 'dryrun'
    dryrun = 1;
  end
 
end

try
  dryrun;
catch
  dryrun = 0;
end

try
  mmc.db = db;
  conn_id = mysql_make_conn(mmc);
catch
  conn_id = mysql_make_conn;
end


%get response table name
sql_resp_table = sprintf(['select response_table from experiment where' ...
		    ' experiment_title = ''%s'''],exp_name);
resp_table_name = mysql(conn_id,sql_resp_table);

resp_table_name = resp_table_name{1};


subjectFields = {'gender','ethnic_category','racial_category'};

nSubFields = length(subjectFields);
for thisSubField = subjectFields

  thisSubField = thisSubField{1};
  
  switch(thisSubField)
    
   case 'gender'
    thisQID = gID;
    
    enum_map = {'F','Female';
		'M','Male';
		'','No Answer'};
    
   case 'ethnic_category'
    thisQID = eID;
    
    enum_map = {'HL','Hispanic or Latino';
		'not_HL','Not Hispanic or Latino';
		'unknown','Unknown (or care not to report)'};
    
   case 'racial_category'
    thisQID = rID;
    enum_map = {'NatAm','American Indian';
		'NatHaw','Native Hawaiian';
		'AfrAm','Black or African American';
		'Asian','Asian';
		'Cauc','Caucasian';
		'More','More than One Race';
		'Unknown','Unknown or No Answer';
		'Other','Other'};
  end
     
  %get responses for this question
  sql_get_resp = sprintf('select subject_id,response_enum from %s where question_id = %d',resp_table_name,thisQID);
  [subID,resp_enum] = mysql(conn_id,sql_get_resp);
  
  %resolve enum values
  sql_get_dfid = sprintf('select answer_format_id from question_x_data_format where question_id = %d',thisQID);
  dfID = mysql(conn_id,sql_get_dfid);
  sql_get_enum_strings = sprintf('select enum_values from data_format where data_format_id = %d',dfID);
  enum_all = mysql(conn_id,sql_get_enum_strings);
  enum_tokens = regexp(enum_all{1},'\"([^,]+)\",?','tokens');
  nenum = length(enum_tokens);

  for ienum = 1:nenum
    enum_strings{ienum} = enum_tokens{ienum}{1};
  end

  resp_vals = enum_strings(log2(resp_enum)+1);
  
  nresp = length(resp_vals);
  for iresp = 1:nresp
    resp_idx = strmatch(resp_vals{iresp},enum_map(:,2));
    mapped_resp_vals{iresp} = enum_map{resp_idx,1};

    %submit to subject table
    sql_update_subject_table = sprintf(['update subject set `%s` = ''%s'' where' ...
		    ' subject_id = ''%s'''],thisSubField, mapped_resp_vals{iresp},subID{iresp});
	
    if(dryrun)
      disp(sql_update_subject_table);
    else
      mysql(conn_id,sql_update_subject_table);
    end
    
  end

  
end
