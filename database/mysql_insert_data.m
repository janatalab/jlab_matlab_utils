function key_list = mysql_insert_data(varargin)
% inserts new records into a database. Returns vector of ids that were inserted
%
% INPUT(S):
% arguments are param/value pairs:
%  'table': name of the table to insert data into
%  'fields': cell array of field names to insert into (should be in the same
%            order as values)
%  'record_values': a cell array of cell arrays. Each inner cell array is a
%                   record with the values to be inserted for a given record, in the same order
%                   specified in 'fields'
%
%  OUTPUT(S):
%  A vector of primary key IDs (tables with composite keys are not supported yet)
%  If a duplicate record was inserted, the highest primary key value is
%  returned, and a warning is given
%
%  Example:
%  inserted_id = mysql_insert_data('table','data_format','fields',{'type','enum_values'},...
%                {{'enum','"yes","no"},{'enum','"left","right"'}});
%
%  Author(s):
%  Sept 23, 2009 - Stefan Tomic, first version
  

for iarg = 1:2:nargin
  
  switch(varargin{iarg})
    
   case 'table'
    table = varargin{iarg+1};
   case 'fields'
    fields = varargin{iarg+1};
   case 'record_values'
    values = varargin{iarg+1};
   case 'conn_id'
    conn_id = varargin{iarg+1};
    
  end
 
end

try
  conn_id;
catch
  conn_id = 0;
end

mysql_check_conn(conn_id,'open');

primary_key = mysql_table_key(table,conn_id);

if(iscell(primary_key))
  error('Composite keys are not supported by this function.');
end

%enclose fieldnames in backquotes
fields = regexprep(fields,'^.*$','`$0`');

fieldString = cell2str(fields,',');

nRecords = length(values);
for iRecord = 1:nRecords
  
  thisRecord = values{iRecord};
    
  nValues = length(thisRecord);
  
  %construct a value string for insert command
  for iValue = 1:nValues

    %only supporting types double or char
    switch(class(thisRecord{iValue}))
     case 'double'
      valueStringArray{iValue} = num2str(thisRecord{iValue});
      whereClause{iValue} = sprintf('%s = %d',fields{iValue},thisRecord{iValue});
     case 'char'
      %note that we're replacing empty strings with NULL. This should be fine,
      %since empty strings are generally represented by NULL in Ensemble
      if(isempty(thisRecord{iValue}))
	valueStringArray{iValue} = 'NULL';
	whereClause{iValue} = sprintf('%s is NULL',fields{iValue});
      else
	valueStringArray{iValue} = ['''' thisRecord{iValue} ''''];
	whereClause{iValue} = sprintf('%s = ''%s''',fields{iValue},thisRecord{iValue});
      end
    end
    
  end

  valueString = cell2str(valueStringArray,',');
  mysql(conn_id,sprintf('insert into %s (%s) values(%s)',table,fieldString,valueString));  

  whereString = cell2str(whereClause,' and ');
  insertedID = mysql(conn_id,sprintf('select %s from %s where %s',primary_key,table,whereString));
  
  %return only one ID per inserted record (multiple IDs may be returned if this
  %is a duplicate record). Throw a warning if this is a duplicate
  if(length(insertedID) > 1)
    warning('Multiple records exist with this set of values');
  elseif(isempty(insertedID))
    error('Failed to insert new record or inserted record not found');
  end
  
  insertedID = insertedID(insertedID == max(insertedID));
  
  key_list(iRecord,1) = insertedID;
  
end



