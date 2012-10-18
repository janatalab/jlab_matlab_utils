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
%  'mysql': mysql info, including conn_id and encryption key if necessary
%
%  OUTPUT(S):
%  A vector of primary key IDs (tables with composite keys are not supported yet)
%  If a duplicate record was inserted, the highest primary key value is
%  returned, and a warning is given
%
%  Example:
%  inserted_id = mysql_insert_data('table','data_format','fields',{'type','enum_values'},...
%                {{'enum','"yes","no"'},{'enum','"left","right"'}});
%
%  Author(s):
%  Sept 23, 2009 - Stefan Tomic, first version
%  03Sep2012 - PJ added support for encrypted fields; removed WHERE clause
%                 handling for returning full matches on inserted IDs.
%                 Matching is performed on primary key
  
encrypted = {};

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
		case 'encrypted_flds'
			encrypted = varargin{iarg+1};
		case 'mysql'
			mysql_params = varargin{iarg+1};
  end
 
end

try
  conn_id;
catch
	try
		conn_id = mysql_params.conn_id;
	catch
		conn_id = 0;
	end
end

% Make sure we have an encryption key if necessary
if ~isempty(encrypted)
	if ~exist('mysql_params','var') || isempty(mysql_params.enc_key)
		error('Encrypted fields requested but no encryption key found')
	end
end

mysql_check_conn(conn_id,'open');

primary_key = mysql_table_key(table,conn_id);

if(iscell(primary_key))
  error('Composite keys are not supported by this function.');
end

%enclose fieldnames in backquotes
origFields = fields;
fields = regexprep(fields,'^.*$','`$0`');
numFields = length(fields);

fieldString = cell2str(fields,',');

% Check to see whether data are coming in as a cell array within a single
% column or a cell array for each of the fields
if size(values,2) == numFields
	dataByCol = true;
	
	% Make sure they are of equal length
	nRecords = cellfun('length',values);
	if any(diff(nRecords))
		error('%s: Columns of values have unequal lengths', mfilename)
	else
		nRecords = nRecords(1);
	end
else
	dataByCol = false;
	nRecords = length(values);
end

key_list = [];
for iRecord = 1:nRecords
	
	if dataByCol
		thisRecord = {};
		for ifld = 1:numFields
			thisRecord{ifld} = values{ifld}{iRecord};
		end
	else
		thisRecord = values{iRecord};
	end
	
	nValues = length(thisRecord);
	
	%construct a value string for insert command
	for iValue = 1:nValues
		encrypt = ismember(origFields{iValue}, encrypted);
		
		if strcmp(origFields{iValue}, primary_key)
			is_primary = true;
			primary_value = thisRecord{iValue};
		else
			is_primary = false;
		end
		
		%only supporting types double or char
		switch(class(thisRecord{iValue}))
			case 'double'
				valueStringArray{iValue} = num2str(thisRecord{iValue});
				if is_primary
					primary_val_str = sprintf('%d', thisRecord{iValue});
				end
				if isempty(thisRecord{iValue})
					valueStringArray{iValue} = 'NULL';
				end
			case 'char'
				%note that we're replacing empty strings with NULL. This should be fine,
				%since empty strings are generally represented by NULL in Ensemble
				if(isempty(thisRecord{iValue}))
					valueStringArray{iValue} = 'NULL';
					whereClause{iValue} = sprintf('%s is NULL',fields{iValue});
				else
					if encrypt
						valueStringArray{iValue} = ['aes_encrypt("' thisRecord{iValue} '","' mysql_params.enc_key '")'];
					else
						valueStringArray{iValue} = ['''' thisRecord{iValue} ''''];
					end
				end
				if is_primary
					primary_val_str = thisRecord{iValue};
				end
		end
		
	end
	
	% Check to see if this particular entry already exists
	mysql_str = sprintf('SELECT %s FROM %s WHERE %s="%s"', primary_key, table, primary_key, primary_val_str);
	existID = mysql(conn_id, mysql_str);
	if ~isempty(existID)
		fprintf('Found %d entries for ID: %s ... skipping ...\n', length(existID), primary_val_str);
		
		% Implement an update protocol if we don't automatically want to skip
		continue
	end
	
	valueString = cell2str(valueStringArray,',');
	sqlString = sprintf('insert into %s (%s) values(%s)',table,fieldString,valueString);
	mysql(conn_id,sqlString);
	
	mysql_str = sprintf('select %s from %s where %s="%s"',primary_key,table,primary_key, primary_val_str);
	insertedID = mysql(conn_id,mysql_str);
	
	%return only one ID per inserted record (multiple IDs may be returned if this
	%is a duplicate record). Throw a warning if this is a duplicate
	if(length(insertedID) > 1)
		warning('Multiple records exist with this set of values');
	elseif(isempty(insertedID))
		error('Failed to insert new record or inserted record not found');
	end
	
	if isnumeric(insertedID)
		insertedID = insertedID(insertedID == max(insertedID));
		key_list(iRecord,1) = insertedID;
	else
		insertedID = insertedID(end);
		key_list{iRecord,1} = insertedID;
	end
	
end



