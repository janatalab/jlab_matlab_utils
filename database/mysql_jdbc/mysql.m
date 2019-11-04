function varargout = jmysql(varargin)
% jmysql.m
%
% Implementation of mysql() via JDBC
%

% 01Nov2019 Petr Janata

% Maintain our driver and connection objects
persistent driver connections

% Get our driver if necessary
if isempty(driver)
    driver = javaObjectEDT('com.mysql.jdbc.Driver');
end

% Initialize our connections array
MAX_CONNECTIONS = 15;
if isempty(connections)
    connections = cell(1,MAX_CONNECTIONS);
end

%
% Our first argument defines our action
%

if nargin < 1
   action = 'status'; 
else
   action = varargin{1};
end

% Check whether our first argument is a connection ID
if ~ischar(action)
    % The action is actually a connection ID
    conn_id = action;
    
    % Our second argument is therefore the action
    action = varargin{2};
    extra_arg_start = 3;
else
    % Default connection ID is zero
    conn_id = 0;
    extra_arg_start = 2;
end

% Get the index of our connection in the array
conn_idx = conn_id+1;

% Determine whether the action is actually a query
possible_actions = {'open','close','status','use','query'};

if ~ismember(action, possible_actions)
    query = action;
    action = 'query';
end

% Execute the appropriate method on our class as a function of the action
switch action
    case 'open'
        % Get the remaining parameters. These have to be sent in a specific
        % order
        host = varargin{extra_arg_start};
        user = varargin{extra_arg_start+1};
        password = varargin{extra_arg_start+2};
        
        % Check whether this connection is already open
        conn = connections{conn_idx};
        if isempty(conn) || conn.isClosed()
            try
                % Set the URL
                url = ['jdbc:mysql://' host '/' ];
                
                % Set the properties
                properties = java.util.Properties();
                properties.setProperty('user', user);
                properties.setProperty('password', password);

                % Create the connection
                conn = driver.connect(url, properties);
                
                % Cache the connection for later use
                connections{conn_idx} = conn;
                
                % Set the output value to be equal to the connection ID
                varargout{1} = conn_id;
                
            catch exception
                throw(MException('jmysql:ConnectionError', char(exception.message)));
            end
        end
        
    case 'close'
        % Fetch the connection from our connections array
        conn = connections{conn_idx};
        
        % Close the connection
        conn.close();
        
        % Remove it from our connections array
        connections{conn_idx} = [];
        
    case 'status'
        conn = connections{conn_idx};
        if ~isempty(conn) && ~conn.isClosed()
            varargout{1} = 0;
        else
            varargout{1} = -1;  % This needs to be changed to reflect actual status codes
        end
        
    case {'use', 'query'}
        if strcmp(action,'use')
            db = varargin{extra_arg_start};
            query = ['use ' db];
        end
       
        % Get our connection
        conn = connections{conn_idx};
        if isempty(conn)
            error('Connection not open')
        end
        
        % Prepare to make a query by opening a statement
        statement = conn.createStatement();
        
        % Execute the query
        if statement.execute(query)
            % Get the result set
            resultset = statement.getResultSet();
        
            % Process the result set
            [outdata, varnames] = process_result_set(resultset);
        
            % Close the result set
            resultset.close()
        else
            outdata = {};
            varnames = {};
        end
        
        % Close the statement
        statement.close();
        
        % Assign the output variables
        if strcmp(action, 'query')
            % In order to be compatible with the C++ version of the mysql
            % connector, we only return the data, not the variable names
            varargout = outdata;
        end
end % switch action

end

function [outdata, varnames] = process_result_set(resultset)
% Get the metaData
metadata = resultset.getMetaData();
ncols = metadata.getColumnCount();

% Count the number of rows
resultset.last(); % scroll to end of data
nrows = resultset.getRow();

% Initialize outputs
varnames = cell(1,ncols);
datatypes = cell(1,ncols);
outdata = cell(1,ncols);

% Get field names and data types
for icol = 1:ncols
   varnames{icol} = char(metadata.getColumnLabel(icol)); 
   datatypes{icol} = char(metadata.getColumnTypeName(icol));
end

% Now get and convert the data
numeric_types = {'DECIMAL','BIGINT UNSIGNED', 'INT UNSIGNED'};
datetypes = {'DATETIME','DATE','TIME'};

% Initialize output data
for icol = 1:ncols
    switch datatypes{icol}
        case numeric_types
            outdata{icol} = nan(nrows,1);
        otherwise
            outdata{icol} = cell(nrows,1);
    end
end

% Scroll back to beginning
if resultset.first()
    have_more_data = true;
    irow = 0;
    
    while have_more_data
        irow = irow+1;
        
        % Loop over fields
        for icol = 1:ncols
            value = resultset.getObject(icol);
            
            if ~resultset.wasNull()
                % Convert the data based on datatype
                switch datatypes{icol}
                    case numeric_types
                        outdata{icol}(irow,1) = double(value);
                    case datetypes
                        outdata{icol}{irow,1} = char(value);
                    otherwise
                        outdata{icol}{irow,1} = value;

                end % switch datatypes{icol}
            end % ~resultset.wasNull()
        end % for icol
        
        % Scroll to next row
        if ~resultset.next()
            have_more_data = false;
        end
    end
end % [outdata, varnames] = process_result_set(resultset)
end

