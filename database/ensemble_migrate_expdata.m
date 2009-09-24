function outData = ensemble_migrate_expdata(varargin)
% migrates experiment data (only questions are currently supported) from one database to another
%
% POSSIBLE FUTURE EXTENSIONS:
%  - migration of complete forms or experiments. This would probably involve implementing recursion
%  - Detection of questions, forms, or experiments that already exist in the
%    destination. In this case, the pre-existing destination IDs in addition to
%    the IDs of records that are created should be returned. This
%    functionality would be essential if implementing the migration of forms
%    and experiments.
%  
% INPUTS
%  the following tag/value pairs:
%        'table' -            the name of the table to migrate (only 'question' table currently supported)
%        'conn_id_from' -     the connection ID of the originating database (optionally can provide 'database_from')
%        'conn_id_to' -       the connection ID of the destination database (optionally can provide 'database_to')
%        'database_from' -    the name of the originating database (used if conn_id_from is not provided)
%        'database_to' -      the name of the destination database (used if conn_id_to is not provided)
%        'primary_key_vals' - a vector (array) of the primary keys corresponding to the records to migrate
%                             (e.g. question_id)
%        'host_from' -        The hostname of the originating database (default host used if ommitted)
%        'host_to' -          The hostname of the destination database (default host used if ommitted)
%        'ignore' -           fields to ignore. They will be ommitted during migration
%                             and the database will assign a default value (if specified)
%
% OUTPUT(S):
%  A vector of primary keys corresponding to the records that were created at
%  the destination.
%
% Example:
%
% destIDs = ensemble_migrate_expdata('table','question','database_from','ensemble_tarp','database_to',...
%                                    'ensemble_main','primary_key_vals',[100 101 102],'ignore','question_category');
%            
%
% **********************************************************************************************************  
%
% "Ensemble" is the proprietary property of The Regents of the University of California ("The Regents.")
%
%  Copyright (c) 2005-09 The Regents of the University of California, Davis campus. All Rights Reserved.
%
%  Redistribution and use in source and binary forms, with or without modification, are permitted by
%  nonprofit, research institutions for research use only, provided the conditions in the included
%  license agreement are met.
%
%  Refer to  for the license agreement,
%  
%
% **********************************************************************************************************
%
% Author(s):
% September 24, 2009 - Stefan Tomic, first version
  
  
for iarg = 1:2:nargin
  
  switch(varargin{iarg})
    
   case 'table'
    table = varargin{iarg+1};
   case 'conn_id_from'
    conn_id_from = varargin{iarg+1};
   case 'conn_id_to'
    conn_id_to = varargin{iarg+1};
   case 'database_from'
    database_from = varargin{iarg+1};
   case 'database_to'
    database_to = varargin{iarg+1};
   case 'primary_key_vals'
    primary_key_vals = varargin{iarg+1};
   case 'host_from'
    host_from = varargin{iarg+1};
   case 'host_to'
    host_to = varargin{iarg+1};
   case 'ignore'
    ignore = varargin{iarg+1};
    
  end
end

% if host_from and host_to were not supplied, then make these blank
% mysql_make_conn will assign the default host
if(~exist('host_from','var'))
  host_from = '';
end

if(~exist('host_to','var'))
  host_to = '';
end

try
  
  %set ignore to only those fields which are supported to ignore
  switch(table)
    case 'question'
     ignore = intersect(ignore,'question_category');
  end
  
catch
  ignore = {};
end

try
  status_conn_id_from = mysql(conn_id_from,'status');
catch
  %if try failed, conn_id variable is probably not set
  status_conn_id_from = 1;

  %find the next available conn id
  for conn_id_from = 1:9
    if(mysql(conn_id_from,'status'))
      break;
    end
    
  end
  
  if(conn_id_from == 9 && ~mysql(conn_id_from,'status'))
    error('Could not find an available database connection');
  end
  
  %in case this was a connection that had an error status (instead of closed
  %status), close it
  mysql(conn_id_from,'close');
  
end

try
  status_conn_id_to = mysql(conn_id_to,'status');
catch
  status_conn_id_to = 1;

  %if try failed, conn_id variable is probably not set
  
  %find the next available conn id
  for conn_id_to = setdiff(1:9,conn_id_from)
    if(mysql(conn_id_to,'status'))
      break;
    end
  end
  
  if(conn_id_to == 9 && ~mysql(conn_id_to,'status'))
    error('Could not find an available database connection');
  end
  
  %in case this was a connection that had an error status (instead of closed
  %status), close it
  mysql(conn_id_to,'close');
  
end

% if status was non-zero, then the connection id is not valid
% (or it wasn't supplied and was assigned above)
% establish a connection here
if(status_conn_id_from)
  conn_id_from = mysql_make_conn(host_from,database_from,conn_id_from);
end

if(status_conn_id_to)
  conn_id_to = mysql_make_conn(host_to,database_to,conn_id_to);
end

switch(table)
 case 'experiment'
  primary_key = 'experiment_id';
  ignoreCompare = {'experiment_id'};
 case 'form'
  primary_key = 'form_id';
 case 'question'
  primary_key = 'question_id';
  ignoreCompare = {'question_id','compqid','data_format_id','locked'};
end

%obtain the metadata
tableMeta = mysql_extract_metadata('table',table,primary_key,primary_key_vals, ...
				 'conn_id',conn_id_from);

%now submit the data to the destination database
%the data types will be different for different tables
%so we need a big switch block in here


  switch(table)
   
   case 'question'
    
    nQuestions = length(primary_key_vals);
    for iQuestion = 1:nQuestions
      oldQid = primary_key_vals(iQuestion);
      metaQids = [tableMeta.question_id];
      metaIdx = find(oldQid == metaQids,1,'first');
      
      %obtaining the question text from the database, since
      %mysql_extract_metadata replaces question_text with the first header
      question_text = mysql(conn_id_from,sprintf('select question_text from question where question_id = %d',oldQid));
      question_text = question_text{1};
      question_category = tableMeta(metaIdx).question_category;
      heading_format = tableMeta(metaIdx).heading_format;
      locked = 'F'; %the question is new to the new database, so unlock it
    
      %mostly, one would just want to ignore question_category, so that's why
      %ignore is only implemented here
      insertFields = setdiff({'question_text','question_category','heading_format','locked'},ignore);
      nflds = length(insertFields);
      for ifld = 1:nflds
	insertVals{1}{ifld}  = eval(insertFields{ifld});
      end
      newQid = mysql_insert_data('conn_id',conn_id_to,'table','question','fields',insertFields,'record_values',insertVals);
  
      qidMap(iQuestion,1:2) = [oldQid newQid];
      
    end
    
    
    
    nMeta = length(tableMeta);
    for iMeta = 1:nMeta
      origQid = tableMeta(iMeta).question_id;
      dfType = tableMeta(iMeta).type;
      enum_values = tableMeta(iMeta).enum_values;
      subquestion = tableMeta(iMeta).subquestion;
      heading = tableMeta(iMeta).heading;
      range = tableMeta(iMeta).range;
      default = tableMeta(iMeta).default;
      html_field_type = tableMeta(iMeta).html_field_type;
      required = tableMeta(iMeta).required;
    
      %try to find a suitable data_format_id in the new database, if this
      %type already exists
    
      %reformat enum_values to values that exist in the database
      if(isempty(enum_values))
	expr_enum_values = 'is NULL';
      else
	%escape any special characters in the enum
	nEnum = length(enum_values);
	for iEnum = 1:nEnum
	  esc_enum_values{iEnum} = regexptranslate('escape',enum_values{iEnum});
	end
	expr_enum_values = regexprep(esc_enum_values,'^.*$','[[.quotation-mark.]]$0[[.quotation-mark.]]');
	expr_enum_values = ['regexp ''^' cell2str(expr_enum_values,'[[:space:]]*,[[:space:]]*') '$'''];
      end
      
      %see if this data type exists in destination
      dfID = mysql(conn_id_to,sprintf('select data_format_id from data_format where type = ''%s'' and enum_values %s',...
				    dfType,expr_enum_values));
    
      if(~isempty(dfID))
	newDfID = dfID;
      else
	if(~isempty(enum_values))
	  quoted_enum_values = regexprep(enum_values,'^.*$','"$0"');
	  enum_str = cell2str(quoted_enum_values,',');
	  newDfID = mysql_insert_data('conn_id',conn_id_to,'table','data_format','fields',{'type','enum_values'},'record_values',{{dfType,enum_str}});
	else
	  newDfID = mysql_insert_data('conn_id',conn_id_to,'table','data_format','fields',{'type'},'record_values',{{dfType}});
	end
      end
   
      destQid = qidMap(find(qidMap(:,1) == origQid),2);
    
      %insert new record to question_x_data_format
      qdfFieldNames = {'question_id','subquestion','answer_format_id','heading','range','default',...
		       'html_field_type','required'};
      qdfValues = {{destQid,subquestion,newDfID,heading,range,default,html_field_type,required}};
      mysql_insert_data('conn_id',conn_id_to,'table','question_x_data_format','fields',qdfFieldNames,'record_values',qdfValues);
      
    end %for iMeta    
    primary_key_list = qidMap(:,2);
    
  end %switch(table)



%if we made the database connections in this function, close them
if(status_conn_id_from)
  mysql(conn_id_from,'close');
end

if(status_conn_id_to)
  mysql(conn_id_to,'close');
end

outData = primary_key_list;
  