%
% Reads the contents of the specified GUI config file
% and puts the results into the return stuct "c" which will become 
% the "config" field of the main GUI data struct. 
% Checks for data validation and setting default values for missing
% variables are all done here.
% 
% To add a new variable, you only need to add it to the config file,
% check for its value here, and assign it to a position within the
% return struct "c". Then the main GUI can do whatever it wants with it. 
%
function [c] = ens_initConfig( fname )

	% if file name doesn't exist then config will come back from parseConfig
	% as an empty list, and all the config values below will be filled with defaults
	config = parseConfig( fname );

	DEFAULT_WIDTH = 750;
	DEFAULT_HEIGHT = 500;
	DEFAULT_BOTTOM_HEIGHT = 48;
	DEFAULT_TOP_HEIGHT = 24;
	DEFAULT_CENTER_BOTTOM_HEIGHT = 26; % height of edit field area under center uitree
	
	DEFAULT_DB_NAME = 'ensemble_main';
	DEFAULT_DB_SERVER = 'atonal';
	DEFAULT_CONN_ID = 9;

	MIN_HEIGHT = 100;
	

	% mysql database settings with default values
	c.db_name = getValue( config, 'db_name' );
	if( ~length( c.db_name ) )
		c.db_name = DEFAULT_DB_NAME;
	end
	c.db_server = getValue( config, 'db_server' );
	if( ~length( c.db_server ) )
		c.db_server = DEFAULT_DB_SERVER;
	end
	c.conn_id = str2double( getValue( config, 'conn_id' ) );
	if( isnan( c.conn_id ) | c.conn_id < 0 )
		c.conn_id = DEFAULT_CONN_ID;
	end


	% GUI positions and dimensions
	c.dim.h = 	str2double( getValue( config, 'height' ) );
	if( isnan(c.dim.h) | c.dim.h <= 0 )
		c.dim.h = DEFAULT_HEIGHT;
	end
	if( c.dim.h < MIN_HEIGHT )
		c.dim.h = MIN_HEIGHT;
	end

	c.dim.w = 	str2double( getValue( config, 'width' ) );
	if( isnan(c.dim.w) | c.dim.w <= 0 )
		c.dim.w = DEFAULT_WIDTH;
	end

	c.dim.bh = 	str2double( getValue( config, 'bottom_height' ) );
	if( isnan(c.dim.bh) | c.dim.bh <= 0 )
		c.dim.bh = DEFAULT_BOTTOM_HEIGHT;
	end

	c.dim.cbh = 	str2double( getValue( config, 'center_bottom_height' ) );
	if( isnan(c.dim.cbh) | c.dim.cbh <= 0 )
		c.dim.cbh = DEFAULT_CENTER_BOTTOM_HEIGHT;
	end

	c.dim.th = 	str2double( getValue( config, 'top_height' ) );
	if( isnan(c.dim.th) | c.dim.th <= 0 )
		c.dim.th = DEFAULT_TOP_HEIGHT;
	end

	% these may change during development but after that
	% these should stay constant
	c.dim.b1h = 20;		% button style 1 height 
	c.dim.b1w = 20;		% button style 1 width 
	
	% size of monitor/resolution
	screenSize = get(0,'ScreenSize');
	c.dim.sh = screenSize(4);
	c.dim.sw = screenSize(3);
	c.dim.minAnTreeHeight = 80;


	checkNum = 1;
	slotNum = 1;
	while(1)
		checkName = sprintf( 'analysisFunction%d', checkNum );
		if( length( getValue( config, checkName ) ) > 0 )
			if( exist( getValue( config, checkName ) ) ~= 2 )
				fprintf( 'Error: analysis function ''%s'' not found in path. Skipping.\n', getValue( config, checkName ) );
				checkNum = checkNum + 1;
				continue;
			end
			c.analysisFunction{slotNum} = getValue( config, checkName );
			checkName = sprintf( 'analysisComment%d', checkNum );
			if( length( getValue( config, checkName ) ) > 0 )
				c.analysisComment{slotNum} = [getValue( config, checkName ) sprintf('\n\nHelp:\n') help( c.analysisFunction{slotNum} )];
			else
				c.analysisComment{slotNum} = [sprintf( 'Help:\n' ) help( c.analysisFunction{slotNum} )];
			end
		else
			break;
		end

		slotNum = slotNum + 1;
		checkNum = checkNum + 1;

	end



function [res] = getValue( config, name )
	res = [];
	
	if( ~length( config ) )
		return;
	end

	for x = 1 : length( config.name ) 

		if( strcmp( name, config.name{x} ) )
			res = config.value{x};
			return;
		end

	end


function [res] = parseConfig( fname )

	res = [];

	f = fopen( fname, 'r' );

	if( f == -1 )
		return;
	end

	x = 1;

	while 1

		l = fgetl(f);

		if( ~ischar(l) )
			break;
		end

		% lines that begin with # are comments
		if( length(l) < 3 | l(1) == '#' )
			continue;
		end

		% locate the | pipe character that separates the
		% config variable name and its value
		p = strfind(l,'|');

		% if there's no |, or if it's the first
		% char (which wouldn't allow a valid name),
		% skip this line. 
		if( length(p) ~= 1 | p < 2 )
			continue;
		end

		name = l(1:p-1);
		value = l(p+1:end);

		res.name{x} = name;
		res.value{x} = value;

		x = x + 1;

	end

	fclose(f);


