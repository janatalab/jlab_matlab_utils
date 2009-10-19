function tableStruct = mysql_extract_metadata(varargin)
% Retrieves experiment, form, question or data_format info.
%
% tableStruct = mysql_extract_metadata(varargin)
%
% examples of usage:
%   tableStruct1 = mysql_extract_metadata('table','experiment','experiment_id',[2 3 4]);
% 
%   tableStruct2 = mysql_extract_metadata('table','experiment','experiment_title',{'Groove' 'autobio_v1'});
%
%   tableStruct3 = mysql_extract_metadata('table','form','form_id',[100 101 102]);
%
%   tableStruct4 = mysql_extract_metadata('table','experiment');  %extracts metadata on all experiments
%
% extracts data from experiment, form, question, or data_format
% tables in the Ensemble database and reorganizes the records into
% arrays of structs, where each field in the struct corresponds to
% a field in the Ensemble table. Tables are traversed in a top-down 
% fashion, where experiments are considered top level and link to
% forms, questions, and data formats.
%
% arguments are tag, value pairs, where possible tags are:
% 
% table: either 'experiment','form','question',or 'data_format'
%
% conn_id: the mysql connection ID to use. If none is specified, 
%          a default conn_id of 1 will be used.
%
% keep_db_open: normally only used for recursively calling this
%               function so that the database isn't closed during
%               successive function calls.
%
% <fieldname>: a field name of the table being extracted. The
%              specified fields do not have to be key fields.
%

% 19 Jan 2007, First Version, S.T.
% 03/15/07 PJ - added conn_id to a couple of mysql calls that were missing
%               it 
% 05/01/07 PJ - enable passing in of host and database arguments for
%               establishing connections to databases other than the
%               default database


extractDataArgs = {};
for(itag = 1:2:length(varargin))
  
  switch(varargin{itag})
    case 'table'
      table = varargin{itag+1};
     
    case 'host'
      host = varargin{itag+1};
      
    case 'database'
      database = varargin{itag+1};
      
    case 'conn_id'
    conn_id = varargin{itag+1};
    
   case 'keep_db_open'
    keep_db_open = varargin{itag+1};
    
   case 'addJunctionTableFields'
    addJunctionTableFields = varargin{itag+1}{:};
    
   case 'addJunctionTableVals'
    addJunctionTableVals = varargin{itag+1}{:};
					   
   otherwise
    extractDataArgs{end+1} = varargin{itag};
    extractDataArgs{end+1} = varargin{itag+1};
    
  end
end


% Check for connection to database
try database(1); 
catch
  database = []; 
end
try host(1); 
catch
  host = []; 
end
try conn_id(1);
catch   
  conn_id = 1;
  disp('making connection to db');
  mysql_make_conn(host,database,conn_id);
  temp_conn_id = 1;
  varargin{end+1} = 'conn_id';
  varargin{end+1} = conn_id;
end

try 
  addJunctionTableFields;
catch
  addJunctionTableFields = {};
end

try 
  addJunctionTableVals;
catch
  addJunctionTableVals = {};
end


switch table
 
 case 'experiment'
  parentKeyField      = 'experiment_id';
  childTable = 'form';
  childKeyField   = 'form_id';
  foreignKeyField    = childKeyField;
  childJunctionTable = 'experiment_x_form';
  orderChildrenBy        = 'form_order';
  junctionTableInfo = {'form_order','form_handler','goto','repeat','condition','condition_matlab','stimulus_matlab','break_loop_button'};
  flattenChildren = 0;
  
 case 'form'
  parentKeyField      = 'form_id';
  childTable = 'question';
  childKeyField   = 'question_id';
  foreignKeyField    = childKeyField;
  childJunctionTable = 'form_x_question';
  orderChildrenBy = 'form_question_num';
  junctionTableInfo = {'question_iteration','form_question_num'};
  flattenChildren = 0;
  
 case 'question'
  parentKeyField      = 'question_id';
  childTable = 'data_format';
  childKeyField   = 'data_format_id';
  foreignKeyField    = 'answer_format_id';
  childJunctionTable = 'question_x_data_format';
  orderChildrenBy = 'subquestion';
  junctionTableInfo = {'subquestion','heading','range','default','html_field_type','required'};
  flattenChildren = 1;
  
 case 'data_format'
  parentKeyField = 'data_format_id';
  childTable = '';
  junctionTableInfo = {};
  flattenChildren = 0;
  
  
 case 'stimulus'
  parentKeyField = 'stimulus_id';
  childTable = 'attribute';
  childKeyField = 'attribute_id';
  foreignKeyField = childKeyField;
  childJunctionTable = 'stimulus_x_attribute';
  orderChildrenBy = 'attribute_id';
  junctionTableInfo = {'attribute_value_double', ...
		    'attribute_value_text'};
  flattenChildren = 0;
  
 case 'attribute'
  parentKeyField = 'attribute_id';
  childTable = '';
  junctionTableInfo = {};
  flattenChildren = 0;
  
 otherwise
  error(sprintf('Unrecognized table: %s',table));
  
end

[tableInfo{1:6}] = mysql(conn_id,sprintf('describe %s',table));
inputFields = extractDataArgs(1:2:end);
[validFields,fieldLocs] = ismember(inputFields,tableInfo{1});
if(any(validFields == 0))
  invalidFieldIdxs = find(validFields == 0);
  for idxInvalidField = invalidFieldIdxs
    disp(sprintf('Invalid Field: %s',inputFields{idxInvalidField}))
    error('Cannot continue since one or more incorrect fields have been entered.');
  end
end


[values,tableFieldNames] = mysql_extract_data('table',table,'conn_id',conn_id,extractDataArgs{:});

%since mysql_extract_data does not retrieve multiple
%records for repeated keys and automatically sorts them, we need to reorder the
%results here

[tf,fldIdx]  = ismember(parentKeyField,tableFieldNames);
[tf,parentKeyLoc] = ismember(parentKeyField,extractDataArgs(1:2:end));

%only reorder if key fields were given. If key fields were not
%given, no reordering is necessary since this would have been done
% at the top level (a selection was made that did not include
% any keys) in which case the function should just return the
% records in the order they were given
if(parentKeyLoc > 0)
  %we found the odd ordered index, now convert to absolute index	
  parentKeyLoc = parentKeyLoc * 2 - 1;
  parentKeyVals = extractDataArgs{parentKeyLoc+1};
  [tf,keyLocs] = ismember(parentKeyVals,values{fldIdx});
  
  if(any(keyLocs == 0))
    
    emptyKeys = find(keyLocs == 0);
    keyLocs(emptyKeys) = [];
  end

  %reorder the values
  for iField = 1:length(tableFieldNames)
    values{iField} = values{iField}(keyLocs);
  end

end

%add junction table tableFieldNames and values if passed through recursion
tableFieldNames = {tableFieldNames{:} addJunctionTableFields{:}};
values     = {values{:}     addJunctionTableVals{:}};

if(strcmp(table,'question'))
  tableFieldNames = {'compqid' tableFieldNames{:}};
  values = {repmat(NaN,size(values{1})) values{:}};
end

%construct arguments for mkstruct
mkstrArgs = cell(1,(length(tableFieldNames)*2));
mkstrArgs(1:2:end) = tableFieldNames;
mkstrArgs(2:2:end) = values;

% make a structure with the tableFieldNames and values
% and convert the structure to an array of structs
tableStruct = convert_structarray(mkstruct(tableFieldNames,mkstrArgs));

%add a field to contain the struct for the child table
if(~isempty(childTable))
  tableStruct(1).(childTable) = [];
end

%initialize newTableStruct (to be used in following loop)
newTableStruct = struct([]);

%go through each record in tableStruct and see if 
%there are children structures that need to be added to this structure
for irecord = 1:length(tableStruct)
  
  if(~isempty(parentKeyField))
    parentKeyVal = getfield(tableStruct(irecord),parentKeyField);
  end

  childStruct = struct([]);
  addTableStructs = struct([]);
  
  %get child information only if we have a value for the parent key
  if(exist('parentKeyVal','var') & ~isempty(parentKeyVal))  
    
 
 
   %childTable contains the name of the table for the child
   %(e.g. 'question' table is the child of the 'form' table
   %obtain child info through the junction table (e.g. form_x_question)
   if(~isempty(childTable))
   
      %construct fieldname strings suitable for sql that do not
      %allow for ambiguity of table
      sqlJunctionTableInfo = strcat([childJunctionTable '.`'],junctionTableInfo(:),'`');
      sqlJunctionTableInfoText = sprintf(',%s', sqlJunctionTableInfo{:});
      sqlForeignKeyFieldText = sprintf('%s.`%s`',childJunctionTable,foreignKeyField);
      sqlChildKeyFieldText = sprintf('%s.`%s`',childTable,childKeyField);
      
      %get the foreign key values to the child table as well as any
      %info fields from the junction table (e.g. question_x_data_format.subquestion)
      sql_junctionTable = sprintf('select %s%s,`%s` from %s join %s on (%s = %s) where %s = %d order by %s',...
				  sqlForeignKeyFieldText,...
				  sqlJunctionTableInfoText,...
				  parentKeyField,...
				  childJunctionTable,...
				  childTable,...
				  sqlChildKeyFieldText,...
				  sqlForeignKeyFieldText,... 
				  parentKeyField,...
				  parentKeyVal,...
				  orderChildrenBy);
     
      %separate the retrieved keys to the child table and the info fields
      junctionTable = cell(1,length(junctionTableInfo(:))+2);
      [junctionTable{:}] = mysql(conn_id,sql_junctionTable);
      childrenKeyVals = junctionTable{1};
      
      if(length(junctionTableInfo) > 0)
	junctionTableInfoVals = junctionTable(2:end-1);
      else
	junctionTableInfoVals = [];
      end
      
      junctionParentKeyVals = junctionTable{end};
      
      %if any child keys were obtained, get metadata on the child tables
      if(length(childrenKeyVals) > 0)
	childStruct = mysql_extract_metadata('table',childTable,...
					     childKeyField,childrenKeyVals,...
					     'conn_id',conn_id,...
					     'keep_db_open',1,...
					     'addJunctionTableFields',{junctionTableInfo},...
					     'addJunctionTableVals',{junctionTableInfoVals});
      end
	
   end %if ~isempty childTable
    
  end %if(~isempty(parentKeyVal)

  if(isempty(addTableStructs))
    addTableStructs = tableStruct(irecord);
  end
 
  
  %if child metadata was not obtained for a question (due to a bad foreign key)
  %get child key fields and construct an
  %empty child struct so that all tableStruct fields match  
  if(strcmp(table,'question') & isempty(childStruct))
    [childInfo{1:6}] = mysql(conn_id,sprintf('describe %s',childTable));
    childFields = childInfo{1};
    childStruct = mkstruct({childFields{:} junctionTableInfo{:}});
  end
  
  %add the child structure if it was obtained
  if(~isempty(childStruct))
    addTableStructs = setfield(tableStruct(irecord),childTable,childStruct);
  end
  
  
  %flatten the children structs with the parent struct if set
  %to true (currently only done for question table)
  if(flattenChildren)
    addTableStructs = flatten_children_structs(addTableStructs);
  end  
  
  
  % add the new Table structs. We need to do this to a variable
  % other than tableStructs (newTableStruct) in order to leave tableStruct intact through
  % each iteration of the loop
  if(isempty(newTableStruct))
    newTableStruct = addTableStructs;
  else
    newTableStruct(end+1:end+length(addTableStructs)) = addTableStructs;
  end
  
end %for

%if a new table structure was created in the for loop above
%then return the new table structure, otherwise return the original
%tableStruct
if(~isempty(newTableStruct))
  tableStruct = newTableStruct;
end

if(strcmp(table,'question'))
  for tsIdx = 1:length(tableStruct)
    
    if(strcmp(table,'question'))
      qid = tableStruct(tsIdx).question_id;
      sqid = tableStruct(tsIdx).subquestion;
      compqid = qid +  sqid ./100;
      tableStruct(tsIdx).compqid = compqid;
      %copy the heading to question_text if the text for the
      %subquestion is located in heading field 
      %(true for subquestion # >=2)
      if(~isempty(tableStruct(tsIdx).heading))
	tableStruct(tsIdx).question_text = tableStruct(tsIdx).heading;
      end
    end
    
  end

end


if (strcmp(table,'data_format'))
  
  for tsIdx = 1:length(tableStruct)
   
    tmpEnumVals = tableStruct(tsIdx).enum_values;
    tokenizedEnums= regexp(tmpEnumVals,'"([^"]*)"','tokens');
    if(~isempty(tokenizedEnums))
      tokenizedEnums = [tokenizedEnums{:}];
      tableStruct(tsIdx).enum_values = {tokenizedEnums{:}};
    end
    
  end
   
end


if(~exist('keep_db_open','var') & exist('temp_conn_id','var'))
  mysql(conn_id,'close');
end

return


