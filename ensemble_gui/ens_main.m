%
% ens_main() is the entry point for the Ensemble analysis management GUI.
% No arguments required. 
% The startup routine looks for a config file called 
%
%
%
%
%
function [mainFig] = ens_main( varargin )
%
% Note: code will look ugly unless viewed in a wide editor window (100-132 chars).
%
% ### we're using "z" as the main guidata variable name for convenience in development.
% this can be changed to something more meaningful but i would actually vote for leaving it as z - JG
%
% Change history:
% 05/07/07: initial handoff of completed v1.0 - JG
%


configFileName = 'ens_config.ini'; % default config file name
arg_conn_id = -1;	%  -1 means no conn_id arg was passed from command line. 
for x = 1 : 2 : length(varargin)
	switch( varargin{x} )
		case 'config'
			configFileName = varargin{x+1};
		case 'conn_id'
			arg_conn_id = varargin{x+1};
	end
end





% loads all values from config file and does some other low level stuff.
% creation of GUI objects is done below. 
% if arg_conn_id is valid, it trumps the default conn_id and/or the conn_id specified in config file
z = init( configFileName, arg_conn_id ); 


if( ~length(z) )
	msg( 'error: init() returned nothing', 'msg' );
	mainFig = [];
	return;
end


% create the main GUI figure
z.handle.window.main = figure(	'tag', 'z.handle.window.main', ...
											'name', 'Ensemble Analysis Job Manager', ...
											'menubar', 'none', ...
											'resize', 'on', ...
											'position', [16 z.config.dim.sh-z.config.dim.h-16 z.config.dim.w z.config.dim.h] );
set( z.handle.window.main,					'resizefcn', {@resizeMainFigureCB, z.handle.window.main} );
mainFig = z.handle.window.main; % store figure handle in function return variable



% create panels
z.handle.panel.left = uipanel(		'units', 'pixels' );
z.handle.panel.midTop = uipanel(		'units', 'pixels' );
z.handle.panel.midCenter = uipanel(	'units', 'pixels' );
z.handle.panel.midBottom = uipanel(	'units', 'pixels' );
z.handle.panel.right = uipanel(		'units', 'pixels' );
z.handle.panel.bottom = uipanel(		'units', 'pixels' );
z.handle.panel.top = uipanel(			'units', 'pixels' );


% create input menus and widgets 


% experiment selection drop-down list for top panel
z.handle.widget.top.selectExperiment = uicontrol( ...
											'style', 'popupmenu', ...
											'string', z.data.meta.expListStr, ...
											'parent', z.handle.panel.top, ...
											'callback', {@selectExperimentListCB, z.handle.window.main} );



% button to add new analysis unit to the chain
z.handle.widget.top.newan = uicontrol( ...
											'style', 'pushbutton', ...
											'parent', z.handle.panel.top, ...
											'string', 'new', ...
											'enable', 'on', ...
											'callback', {@newAnalysisCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.top.newan' );


% button to delete selected job unit from the list
z.handle.widget.top.deletean = uicontrol( ...
											'style', 'pushbutton', ...
											'parent', z.handle.panel.top, ...
											'string', 'delete', ...
											'enable', 'off', ...
											'callback', {@deleteAnalysisCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.top.deletean' );


% button to quit
z.handle.widget.top.quit = uicontrol( ...
											'style', 'pushbutton', ...
											'parent', z.handle.panel.top, ...
											'string', 'quit', ...
											'enable', 'on', ...
											'callback', {@quitCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.top.quit' );


% buttons to reposition analysis units in the list
z.handle.widget.bottom.movedown = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', '+', ...
											'enable', 'off', ...
											'callback', {@moveCB, z.handle.window.main}, ... 
											'tag', 'movedown' );
z.handle.widget.bottom.moveup = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', '-', ...
											'enable', 'off', ...
											'callback', {@moveCB, z.handle.window.main}, ... 
											'tag', 'moveup' );


% button to add left panel selections to filt
z.handle.widget.bottom.addinclude = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', 'include', ...
											'enable', 'off', ...
											'callback', {@addToFilterCB, z.handle.window.main}, ... 
											'tag', 'addinclude' );

% button to subtract left panel selections to filt
z.handle.widget.bottom.addexclude = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', 'exclude', ...
											'enable', 'off', ...
											'callback', {@addToFilterCB, z.handle.window.main}, ... 
											'tag', 'addexclude' );


% save job spec
z.handle.widget.bottom.save = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.top, ...
											'string', 'save', ...
											'enable', 'off', ...
											'callback', {@saveCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.top.save' );
% save with results attached?
z.handle.widget.bottom.includeResults = uicontrol( ...	
											'style', 'togglebutton', ...
											'parent', z.handle.panel.top, ...
											'enable', 'off', ...
											'tag', 'z.handle.widget.top.includeResults' );

% load job spec
z.handle.widget.bottom.load = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.top, ...
											'string', 'load', ...
											'enable', 'on', ...
											'callback', {@loadCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.top.load' );


% execute the job spec 
z.handle.widget.bottom.go = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', '!', ...
											'enable', 'off', ...
											'callback', {@goCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.bottom.go' );


% this button will dump debug info to the std matlab output
z.handle.widget.bottom.dump = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', '?', ...
											'enable', 'on', ...
											'visible', 'off', ...
											'callback', {@dumpCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.bottom.dump' );
if( 1 ) % DEBUG
	set( z.handle.widget.bottom.dump, 'visible', 'on' );
end


% this button will display info about the selected analysis unit
z.handle.widget.bottom.anInfo = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', 'info', ...
											'enable', 'off', ...
											'callback', {@anInfoCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.bottom.anInfo' );

% this button will set the analysis of the currently selected unit to the
% currently selected analysis
z.handle.widget.bottom.anAdd = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.bottom, ...
											'string', 'apply', ...
											'enable', 'off', ...
											'callback', {@anApplyCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.bottom.anAdd' );


% lists all analysis units selected for the current spec
z.handle.widget.right.analysisListbox = uicontrol( ...
											'style', 'listbox', ...
											'parent', z.handle.panel.right, ...
											'string', '', ...
											'value', 0, ...
											'callback', {@analysisSelectCB, z.handle.window.main}, ...
											'tag', 'z.handle.widget.right.analysisListbox' );

% populate analysis list for right-window
% the analysis info has already been extracted from the database in the init() routine.
% this list is static for the life of this GUI so we only have to set it up once here. 
if( length( z.data.meta.anList{2} ) )
	s = cell2str(z.data.meta.anList{2}(1));
	for x = 2 : length( z.data.meta.anList{2} )
		s = [s '|' cell2str(z.data.meta.anList{2}(x))];
	end
	set( z.handle.widget.right.analysisListbox, 'value', 1 );
	set( z.handle.widget.right.analysisListbox, 'string', s );
	set( z.handle.widget.bottom.anInfo, 'enable', 'on' );
end


% lists all analysis units. initially empty
z.handle.widget.mid.jobSpecListbox = uicontrol( ...
											'style', 'listbox', ...
											'parent', z.handle.panel.midTop, ...
											'string', '', ...
											'value', 0, ...
											'callback', {@regenerateTreeForSelectedJob, z.handle.window.main}, ...
											'tag', 'z.handle.widget.mid.jobSpecListbox' );


% create all the input widgets for the analysis unit fields appearing in the center-bottom section.
% these will all be invisible initially and made visible as necessary for editing- only
% 1 will ever be enabled and visible at any given time, since they all will have the same position.
z.handle.widget.mid.fieldEditString = uicontrol( ...
											'style', 'edit', ...
											'string', '', ...
											'visible', 'off', ...
											'enable', 'off', ...
											'callback', {@fieldEditTextCB, z.handle.window.main}, ...
											'parent', z.handle.panel.midBottom );
z.handle.widget.mid.fieldEditPopupmenu = uicontrol( ...
											'style', 'popupmenu', ...
											'string', '<none>', ...
											'visible', 'off', ...
											'enable', 'off', ...
											'callback', {@fieldEditPopupCB, z.handle.window.main}, ...
											'parent', z.handle.panel.midBottom );
z.handle.widget.mid.filtDeleteButton = uicontrol( ...	
											'style', 'pushbutton', ...
											'parent', z.handle.panel.midBottom, ...
											'string', 'delete', ...
											'enable', 'off', ...
											'visible', 'off', ...
											'callback', {@filtDeleteCB, z.handle.window.main}, ... 
											'tag', 'z.handle.widget.mid.filtDeleteButton' );



% these are all plain read-only text fields
z.handle.text.reorder = uicontrol ( ...
											'style', 'text', ...
											'string', 'reorder', ...
											'parent', z.handle.panel.bottom );

z.handle.text.includeResults = uicontrol ( ...
											'style', 'text', ...
											'string', '<-include results', ...
											'parent', z.handle.panel.top );

% status message at bottom of GUI
z.handle.text.status = uicontrol ( ...
											'style', 'text', ...
											'string', 'Welcome to Ensemble.', ...
											'parent', z.handle.panel.bottom );

% variable name in center bottom panel (for editing analysis unit fields)
z.handle.text.fieldName = uicontrol ( ...
											'style', 'text', ...
											'string', '<no edit>', ...
											'parent', z.handle.panel.midBottom );




% create certain nested structs/fields here that wouldn't otherwise
% exist in the beginning, because: 
% there are certain sections of the code that refer to these fields.
% if they refer to the fields before they exist, matlab
% will break. those sections could check for the existence of these
% fields every time it refers to them with isfield() or whatever, but it's
% so much easier and cleaner to just guarantee
% that these fields will always exist, even though they may be == [].
z.data.left.temp = 1;
z.data.left = rmfield( z.data.left, 'temp' ); 
z.data.an{1} = 1;
z.data.an(1) = [];
z.data.left.selected = [];


% store the master struct as guidata for main figure
guidata( z.handle.window.main, z );



% this will always be called last when changes are made
% which affect the layout or size of things - including window resizing.
positionGUIElements(z);






% =========================================
function [] = resizeMainFigureCB( h, v, mainFig )
	z = guidata(mainFig);
	newp = get(z.handle.window.main,'position');
	z.config.dim.w = newp(3);
	z.config.dim.h = newp(4);
	positionGUIElements(z);
	guidata(z.handle.window.main,z); 






% ==================================
% anything that is affected by resizing of the main figure
% needs to be dealt with here. this also initializes the positions at startup.
function [] = positionGUIElements(z)

	ph = z.config.dim.h-z.config.dim.bh-z.config.dim.th-2; % panel height
	otl = floor( z.config.dim.w/3 ); % one thirds left
	ttl = 2 * floor( z.config.dim.w/3 ); % two thirds left or talk to you later

	h = z.handle;

	% find out height of the top middle panel which contains the analysis listbox
	% if there are no units at all, pretend listSize is 1 so the box is visible
	listSize = max([1 length( z.data.an )]);
	
	% set minimum height for the analysis item display in center mid panel.
	% this is specified in the init() fucntion, but don't let it cover 	
	% the entire panel if it's been resized to a smaller height than
	% the config dimension says.
	minAnHeight = z.config.dim.minAnTreeHeight;
	if( minAnHeight > ph-(z.config.dim.cbh+40) )
		minAnHeight = ph-(z.config.dim.cbh+40);
	end

	listHeight = listSize * 16;
	listHeight = min([listHeight ph-(minAnHeight+z.config.dim.cbh)]);
	listHeight = max([listHeight 20]);
	% ------------------------------------------------------------


	if( isfield(h,'panel') )
		set( z.handle.panel.left,	'position', [1 z.config.dim.bh+1 otl-1 ph] );
		set( z.handle.panel.midTop,'position', [otl+1 z.config.dim.bh+(ph-(listHeight+7)) otl-1 listHeight+8] );
		set( z.handle.panel.midCenter,'position',[otl+1 z.config.dim.bh+z.config.dim.cbh+1 otl-1 (ph-listHeight)-z.config.dim.cbh-8] );
		set( z.handle.panel.midBottom,'position',[otl+1 z.config.dim.bh+1 otl-1 z.config.dim.cbh] );
		set( z.handle.panel.right,	'position', [ttl+1 z.config.dim.bh+1 otl-1 ph] );
		set( z.handle.panel.bottom,'position', [1 1 z.config.dim.w-1 z.config.dim.bh] );
		set( z.handle.panel.top,	'position', [1 z.config.dim.bh+ph+2 z.config.dim.w-1 z.config.dim.th] );
	end


	if( isfield( h, 'text' ) )
		set( z.handle.text.status,	'position', [1 1 z.config.dim.w-4 17] );
		set( z.handle.text.fieldName,	'position', [1 3 floor(otl/3)-2 17] );
		set( z.handle.text.includeResults, 'position', [(2*otl)+23 1 116 17] );
		set( z.handle.text.reorder, 'position', [otl+42 z.config.dim.bh-z.config.dim.b1h-6 50 17] );
	end


	if( isfield( h, 'widget' ) )
		set( z.handle.widget.bottom.movedown, 'position', [otl z.config.dim.bh-z.config.dim.b1h-5 20 20] );
		set( z.handle.widget.bottom.moveup, 'position', [otl+21 z.config.dim.bh-z.config.dim.b1h-5 20 20] );
		set( z.handle.widget.bottom.addinclude, 'position', [floor(otl/2)-60 ...
												z.config.dim.bh-z.config.dim.b1h-5 60 z.config.dim.b1h] );
		set( z.handle.widget.bottom.addexclude, 'position', [floor(otl/2)+1 ...
												z.config.dim.bh-z.config.dim.b1h-5 60 z.config.dim.b1h] );
		set( z.handle.widget.bottom.save, 'position', [(2*otl)-40 1 40 20] );
		set( z.handle.widget.bottom.includeResults, 'position', [(2*otl)+2 1 20 20] );
		set( z.handle.widget.bottom.load, 'position', [otl+2 1 40 20] );
		set( z.handle.widget.bottom.go, 'position', [(3*otl/2)-15 z.config.dim.bh-z.config.dim.b1h-5 30 20] );
		set( z.handle.widget.bottom.dump, 'position', [z.config.dim.w-24 z.config.dim.bh-22 20 18] );
		set( z.handle.widget.bottom.anInfo, 'position', [z.config.dim.w-(otl/2)-40+46 z.config.dim.bh-22 34 18] );
		set( z.handle.widget.bottom.anAdd, 'position', [z.config.dim.w-(otl/2)-41 z.config.dim.bh-22 46 18] );

		set( z.handle.widget.top.selectExperiment, 'position', [1 1 otl-4 20] );
		set( z.handle.widget.top.newan, 'position', [(1.5*otl)+5 1 35 20] );
		set( z.handle.widget.top.deletean, 'position', [(1.5*otl)-45 1 50 20] );
		set( z.handle.widget.top.quit, 'position', [z.config.dim.w-37 1 35 20] );

		set( z.handle.widget.mid.jobSpecListbox, 'position', [1 2 otl-5 listHeight] );

		set( z.handle.widget.right.analysisListbox, 'position', [1 1 otl-4 ph-4] );

		set( z.handle.widget.mid.fieldEditString, 'position', [floor(otl/3) 3 floor((2*otl)/3)-4 17] ); 
		set( z.handle.widget.mid.fieldEditPopupmenu, 'position', [floor(otl/3) 3 floor((2*otl)/3)-4 20] ); 
		set( z.handle.widget.mid.filtDeleteButton, 'position', [floor(otl/2)+10 1 50 20] ); 

	end




% =====================================================================
% ### do we want to add functionality to save window position and size?
%
function [] = loadCB( a1, a2, fig )

	z = guidata( fig );

	pos = get( fig, 'position' );
	x = pos(1) + floor(z.config.dim.w/3);
	y = z.config.dim.sh - pos(2) - (z.config.dim.h - z.config.dim.th + z.config.dim.bh);

	[fname, pname] = uigetfile( '*.job', 'location', [x y] );
	fullname = fullfile( pname, fname );

	% matlab won't let you select a file that is a directory or doesn't exist.
	% so all we have to worry about is that it has .job as the extension
	if( length( fname ) < 5 )
		msg( 'no load file specified', 'status' );
		return;
	end		
	if( length( strfind( fname, '.job' )) ~= 1  )
		msg( 'wrong file type  error: GH54', 'errbox' ); 
		msg( sprintf( 'the file you specified is not a .job file: %s', fullname ), 'status' ); 
		return;
	end

	% -mat forces matlab to treat the .job file type as if it were a .mat file
	% old will be a struct with one element: z
	% the data that was present when the file was saved will only
	% be present if the checked the box to include results in the save.
	% old.z.an{x}.results  could potentially be very large in that case. 
	old = load( '-mat', fullname );

	if( isfield( old, 'z' ) )

		data = old.z.data;
		if( isfield( data, 'an' ) )
			z.data.an = old.z.data.an;
			% need to do this here, before it's rendered, in case the list is currently empty.
			% if it is, then jobSpecListBox's value is 0 and matlab will complain 
			set( z.handle.widget.mid.jobSpecListbox, 'value', 1 );

			% if a previous middle tree exists, delete it and its container --------
			% NOTE: this code is copied from jobSpecSelectCB- to make it a function,
			% we'd have to save and load guidata() here again to pass the info back
			% and forth which is costly so as long as it's only duplicated
			% once it's ok for now. 
			h = z.handle;
			if( isfield( h, 'container' ) )
				c = h.container;
				if( isfield( c, 'midTree' ) )
					delete( z.handle.container.midTree );
					z.handle.container = rmfield( z.handle.container, 'midTree' );
				end
			end
			% ----------------------------------------------------------------------


			% ### might be nice to automatically generate the uitree for the 1st unit in the loaded list

		else
			msg( 'your file loaded, but no analysis units were found within', 'msgbox' );	
		end

	else
			msg( 'your file loaded, but no z data was found', 'msgbox' );	
	end


	% if we have nothing in the job list for whatever the reason, disable some buttons
	% but since we can't save an empty analysis list, we should
	% never be able to load one, so length(a.data.an) should always be > 0 here
	if( length( z.data.an ) )
		set( z.handle.widget.bottom.save, 'enable', 'on' );
		set( z.handle.widget.bottom.includeResults, 'enable', 'on' );
		set( z.handle.widget.top.deletean, 'enable', 'on' );
		set( z.handle.widget.bottom.anAdd, 'enable', 'on' );
		
		% if all analysis units have callback functions specified, enable execute button
		set( z.handle.widget.bottom.go, 'enable', 'on' );
		for x = 1 : length( z.data.an )
			if( ~length( z.data.an{x}.fun ) )
				set( z.handle.widget.bottom.go, 'enable', 'off' );
			end
		end
		
	else
		set( z.handle.widget.bottom.save, 'enable', 'off' );
		set( z.handle.widget.bottom.includeResults, 'enable', 'off' );
		set( z.handle.widget.top.deletean, 'enable', 'off' );
		set( z.handle.widget.bottom.go, 'enable', 'off' );
		set( z.handle.widget.bottom.anAdd, 'enable', 'off' );
	end

	% allow reordering position within unit list if there's at least 2
	if( length( z.data.an ) > 1 )
		set( z.handle.widget.bottom.moveup, 'enable', 'on' );
		set( z.handle.widget.bottom.movedown, 'enable', 'on' );
	else
		set( z.handle.widget.bottom.moveup, 'enable', 'off' );
		set( z.handle.widget.bottom.movedown, 'enable', 'off' );
	end
		

	msg( sprintf( 'loaded save file:  %s\n', fullfile( pname, fname) ), 'status' );
	guidata( fig, z );

	updateUnitListDisplay( fig );

	positionGUIElements( z );





function [] = jg( i )

	fprintf( 'jg here %d\n', i );




% ==========================================
function [] = saveCB( a1, a2, fig )

	z = guidata( fig );

	pos = get( fig, 'position' );
	x = pos(1) + floor(z.config.dim.w/3);
	y = z.config.dim.sh - pos(2) - (z.config.dim.h - z.config.dim.th + z.config.dim.bh); % + z.config.dim.bh + 10


	[fname, pname] = uiputfile( '.job', 'save current state', 'location', [x y] );
	fullname = fullfile( pname, fname );

	% this happens when user presses cancel button
	if( length( fname ) < 2 | length( pname ) < 2 )
		msg( 'save cancelled', 'status' );
		return;
	end

	% otherwise, check for invalid input
	if( length( fname ) < 5 | length( strfind( fname, '.job' )) ~= 1 )
		msg( sprintf( 'couldn''t save file: %s', fullname ), 'status' ); 
		return;
	else

		% ### why not just save z.data.an directly since that's all we're keeping?
		% originally there was other stuff in z we needed, and maybe in the future there
		% will be some other settings we want to save like window size/position. 
		anInfo = z.data.an;
		
		% check the includeResults checkbox - if not checked, make sure results are all []
		if( get( z.handle.widget.bottom.includeResults, 'value' ) == 0 )
			for x = 1 : length( anInfo )
				anInfo{x}.results = {}; % ### empty cell array or empty struct?
			end
			% fprintf( 'deleted all existing results from each analysis unit before saving\n' );
		else
			% fprintf( 'including results with save structures\n' );
		end
	

		clear z;
		z.data.an = anInfo;
		save( fullfile( pname, fname ), 'z' );
		msg( sprintf( 'saved: %s', fullfile( pname, fname ) ), 'status' ); 
	end



% ====================================
% this init() is mostly to read params from the config file and load experiment
% meta data from the ensemble mysql database specified in the config file. 
% creation of ui objects, panels, and widgets etc. is not done here. 
function [z] = init( configFileName, arg_conn_id )

	if( exist( configFileName ) ~= 2 )	
		msg( sprintf( '\n\n** Config file "%s" doesn''t exist. using defaults **\n\n', configFileName ) );
	end

	z.config = ens_initConfig( configFileName );

	% if any of these are true, then initConfig had trouble
	if( ~length(z.config) | ~isfield( z.config, 'db_server' ) | ~isfield( z.config, 'db_name' ) | ~isfield( z.config, 'conn_id' ) )
		msg( 'initConfig had problems. err: WD23', 'errbox' );
		z = [];
		return;
	end
		
	
	% if conn_id was passed in at the command line, it trumps the conn_id resulting from initConfig()
	if( arg_conn_id > -1 )
		z.config.conn_id = arg_conn_id;
	end

	mysql_make_conn( z.config.db_server, z.config.db_name, z.config.conn_id ); 


	% load master list of experiments in this database
	z.data.meta.expList = mysql_extract_data(	'table', 'experiment', ...
															'conn_id',  z.config.conn_id, ...
															'extract_flds', {'experiment_id','experiment_title'} );
	if( length( z.data.meta.expList{1} ) )
		z.data.meta.expListStr = cell2str(z.data.meta.expList{2}(1));
		for x = 2 : length( z.data.meta.expList{2} )
			z.data.meta.expListStr = [z.data.meta.expListStr '|' cell2str(z.data.meta.expList{2}(x))];
		end
	else
		z.data.meta.expListStr = 'no experiments available';
	end


	z.data.meta.anList = mysql_extract_data(	'table', 'analysis', ...
															'conn_id',  z.config.conn_id, ...
															'extract_flds', {'analysis_id','function_name','comment'} );
	
	if( isfield( z.config, 'analysisFunction' ) == 1 )
		l = length( z.data.meta.anList{2} );
		for i = 1 : length( z.config.analysisFunction )
			
			z.data.meta.anList{1}(l+i) = l+i; % i don't think this field really matters
			z.data.meta.anList{2}(l+i) = z.config.analysisFunction(i);
		
			% comments should always be filled in with "n/a" or whatever even if
			% the user doesn't specify a comment in the ini file. ens_initConfig
			% should make sure of that. but check here anyway
			if( isfield(  z.config, 'analysisComment' ) == 1 && length( z.config.analysisComment ) >= i )
				z.data.meta.anList{3}(l+i) = z.config.analysisComment(i);
			else
				z.data.meta.anList{3}(l+i) = 'no comment';
			end
		end
	end




% =======================================
function [] = anApplyCB( inarg1, inarg2, fig )
	z = guidata(fig);

	% first make sure we have an anlysis selected on the right,
	% and a created unit selected in the middle. 
	
	rightValue = get( z.handle.widget.right.analysisListbox, 'value' );
	midValue = get( z.handle.widget.mid.jobSpecListbox, 'value' );

	if( midValue < 1 | midValue > size(z.data.an) )
		msg( 'invalid analysis unit selection (middle pane)', 'status' );
		return;
	end

	if( rightValue < 1 | rightValue > length( z.data.meta.anList{2} ) )
		msg( 'invalid analysis unit selection (right pane)', 'status' );
		return;
	end

	try
		z.data.an{midValue}.fun = str2func( cell2str( z.data.meta.anList{2}(rightValue) ) );
	catch
		msg( sprintf( 'error QA32: trying to assign function handle: %s',z.data.meta.anList{2}(rightValue)), 'errbox' );
		return;
	end

	% in case the returned params struct doesn't come with a default filter,
	% save the current filter to reattach to the new params struct. 
	saveFilt = z.data.an{midValue}.params.filt;

	% analysis functions are required to support  fun( 'getDefaultParams' ) 
	evalStr = sprintf( 'z.data.an{midValue}.params = %s( ''getDefaultParams'' );', func2str( z.data.an{midValue}.fun ) );
	try
		eval( evalStr );
	catch
		msg( sprintf( 'error when calling %s( ''getDefaultParams'' )\nsetting params to []', ...
																								func2str( z.data.an{midValue}.fun ) ), 'msgbox' );
		% JG should maybe separate this into its own getDefaultParams function,
		% for cases where the analysis function doesn't supply its own getDefaultParams?

	end

	if( isfield( z.data.an{midValue}.params, 'filt' ) )
		% do nothing- keep the default filter returned by the function
	else
		% populate the filt field with whatever was in there before- might be [] or whatever. 
		z.data.an{midValue}.params.filt = saveFilt;
	end

	% if all analysis units in the list have callback functions specified, enable the execute button
	% assume they all do
	set( z.handle.widget.bottom.go, 'enable', 'on' );
	for x = 1 : length( z.data.an )
		if( ~length( z.data.an{x}.fun ) )
			% if at least one unit doesn't have a CB function, disable execute button
			set( z.handle.widget.bottom.go, 'enable', 'off' );
		end
	end		


	guidata( fig, z );
	
	% this will force the middle pane to reload the uitree for the selected anUnit
	regenerateTreeForSelectedJob( '', '', fig )

	msg( sprintf( 'function for selected analysis unit set to right panel list item %d', rightValue ), 'status' );




% =======================================
function [] = anInfoCB( inarg1, inarg2, fig )
	z = guidata(fig);

	v = get( z.handle.widget.right.analysisListbox, 'value' );
	msg( cell2str( z.data.meta.anList{3}(v) ), 'msgbox' );
 




function [] = quitCB( inarg1, inarg2, fig )

delete( fig );






% ===========================================
function [] = moveCB( inarg1, inarg2, fig )
	z = guidata( fig );


	if( length( z.data.an ) < 2 )
		msg( 'nothing to move', 'status' );
		return;
	end

	tag = get( gcbo, 'tag' );
	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	lim = length( z.data.an );

	if( strcmp( tag, 'movedown' ) )
		if( whichAnUnit < lim & whichAnUnit > 0 )
			temp = z.data.an{whichAnUnit+1};
			z.data.an{whichAnUnit+1} = z.data.an{whichAnUnit};
			z.data.an{whichAnUnit} = temp;
			set( z.handle.widget.mid.jobSpecListbox, 'value', whichAnUnit+1 );
		end
	end

	if( strcmp( tag, 'moveup' ) )
		if( whichAnUnit <= lim & whichAnUnit > 1 )
			temp = z.data.an{whichAnUnit-1};
			z.data.an{whichAnUnit-1} = z.data.an{whichAnUnit};
			z.data.an{whichAnUnit} = temp;
			set( z.handle.widget.mid.jobSpecListbox, 'value', whichAnUnit-1 );
		end
	end

	guidata( fig, z ) 

	% this is not really necessary- we could just call updateUnitListDisplay to see the new list in the
	% new order. but if there's a ui tree open for the current unit, and the requires field
	% is selected, there will be a popup list already created which would have the old
	% list ordering. if you were to then try to modify the requires field by choosing
	% a different unit from the stale popup list, the results are unpredictable. 
	regenerateTreeForSelectedJob( '', '', fig );




% ===========================================
function [] = goCB( inarg1, inarg2, fig )
	z = guidata( fig );

	% make sure each unit has a function that is visible to matlab.
	for x = 1 : length( z.data.an )
		if( exist( func2str( z.data.an{x}.fun ) ) ~= 2 )
			msg( sprintf( 'Function %s for analysis unit %d is not visible to matlab', ...
																func2str( z.data.an{x}.fun ), x ), 'msgbox' );	
			return;	
		end
	end

	params.ensemble.conn_id = z.config.conn_id;
	params.ensemble.db_name = z.config.db_name;
	params.ensemble.db_server = z.config.db_server;
	params.run_analyses = [1:length(z.data.an)];

	msg( sprintf( 'executing analyses%s', sprintf( ' %d', params.run_analyses ) ), 'status' );
	z.data.an = ensemble_jobman( z.data.an, params );

	% save results
	guidata( fig, z );

	regenerateTreeForSelectedJob( '', '', fig );




% ===========================================
function [] = dumpCB( inarg1, inarg2, fig )
	z = guidata( fig );


	msg( 'status dumped to matlab output', 'status' );
	fprintf( '\n\n---- status ------------------------------------\n' );
	fprintf( 'z\n' );

	% it probably isn't necessary to check for these basic stucts
	% every time, but there are cases where buttons and callbacks
	% are executed before some of them are filled in, causing an error.
	% we could keep track of which ones are always explicitly set
	% in the startup routine, but it's hard to keep it all straight. 
	d = z.data;
	fprintf( '\tdata\n' );

	if( length( z.data.an ) < 1 )
		fprintf( '\t\tan = []\n' );
	else
		for x = 1 : length( z.data.an )
			fprintf( '\t\tan{%d}\n', x );
			fprintf( '\t\t\tname = %s\n', z.data.an{x}.name );
			%fprintf( '\t\t\trequires = %s\n', z.data.an{x}.requires{1}.name );
		end
	end

	fprintf( '\t\tleft\n' ); 
	l = z.data.left;
	if( isfield( l, 'selected' ) )
		fprintf( '\t\t\tselected: %d items\n', size( z.data.left.selected,1) );
	else
		fprintf( '\t\t\tselected = null\n' );
	end


	% these are all guaranteed to exist because we create these widgets at startup
	fprintf( '\thandle\n' );
	fprintf( '\t\twidget\n' );
	fprintf( '\t\t\tmid\n' );
	fprintf( '\t\t\t\tjobSpecListbox\n' );
	fprintf( '\t\t\t\t\tvalue = %d\n', get( z.handle.widget.mid.jobSpecListbox, 'value' ) );
	fprintf( '\t\t\tright\n' );
	fprintf( '\t\t\t\tanalysisListbox\n' );
	fprintf( '\t\t\t\t\tvalue = %d\n', get( z.handle.widget.right.analysisListbox, 'value' ) );
	

	fprintf( '\n\n' );







% =========================================
function [] = analysisSelectCB( arg1, arg2, fig )

	% nothing to do here


	


% ==============================
function [] = regenerateTreeForSelectedJob( arg1, arg2, fig )

	z = guidata( fig );

	if( length( z.data.an ) < 1 )
		return;
	end

	whichSelected = get( z.handle.widget.mid.jobSpecListbox, 'value' );

	if( whichSelected < 1 | whichSelected > length( z.data.an ) )
		msg( sprintf( 'selected analysis unit %d is out of bounds', whichSelected ), 'errbox' );
		return;
	end

	% if a previous middle tree exists, delete it and its container --------
	h = z.handle;
	if( isfield( h, 'container' ) )
		c = h.container;
		if( isfield( c, 'midTree' ) )
			delete( z.handle.container.midTree );
			z.handle.container = rmfield( z.handle.container, 'midTree' );
		end
	end
	% ----------------------------------------------------------------------

	z.handle.widget.mid.anTree = getMiddleTree( z.data.an{whichSelected}, z.handle.window.main );
	z.handle.container.midTree = get(z.handle.widget.mid.anTree, 'uicontainer');
	set(z.handle.container.midTree, 'parent', z.handle.panel.midCenter );


	% after diabling any previous controls, enable the textfield for editing name directly even though
	% the uitree hasn't been clicked on. this is sort of a luxury feature
	disableAllfieldEditWidgets( fig );
	z.data.center.vName = 'name';
	z.data.center.vValue = z.data.an{whichSelected}.name;
	set( z.handle.widget.mid.fieldEditString, 'visible', 'on' ); 
	set( z.handle.widget.mid.fieldEditString, 'enable', 'on' ); 
	set( z.handle.widget.mid.fieldEditString, 'string', z.data.center.vValue );
	set( z.handle.text.fieldName, 'string', z.data.center.vName );
	

	guidata( fig, z ) 
	updateUnitListDisplay( fig );
	positionGUIElements( z );





% ======================================
function [] = deleteAnalysisCB( arg1, arg2, fig )

	z = guidata( fig );

	anVal = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	if( anVal < 1 | anVal > length( z.data.an ) )
		return;
	end

	% first delete the display tree for this analysis unit
	h = z.handle;
	if( isfield( h, 'container' ) )
		c = h.container;
		if( isfield( c, 'midTree' ) )
			delete( z.handle.container.midTree );
			z.handle.container = rmfield( z.handle.container, 'midTree' );
			
		end
	end

	% we need this below
	oldName = z.data.an{anVal}.name;

	% now remove the analysis struct from the cell array.
	z.data.an(anVal) = [];

	% disallow reordering of position within unit list if there's not at least 2 of them
	if( length( z.data.an ) < 2 )
		set( z.handle.widget.bottom.moveup, 'enable', 'off' );
		set( z.handle.widget.bottom.movedown, 'enable', 'off' );
	end

	% set the value to 1 by default or 0 if there are no units left. 
	% disable the delete button if there's nothing left to delete.
	if( length( z.data.an ) < 1 )
		set( z.handle.widget.top.deletean, 'enable', 'off' );
		set( z.handle.widget.bottom.save, 'enable', 'off' );
		set( z.handle.widget.bottom.includeResults, 'enable', 'off' );
		set( z.handle.widget.mid.jobSpecListbox, 'value', 0 );
		set( z.handle.widget.bottom.go, 'enable', 'off' );
		set( z.handle.widget.bottom.addexclude, 'enable', 'off' );
		set( z.handle.widget.bottom.addinclude, 'enable', 'off' );
		set( z.handle.widget.bottom.anAdd, 'enable', 'off' );
	else
		% default selection to item 1 after each deletion
		set( z.handle.widget.mid.jobSpecListbox, 'value', 1 );

		% if one of the other units required the deleted unit, 
		% throw up a message and erase the stale dependency
		for x = 1 : length( z.data.an )
			for y = 1 : length( z.data.an{x}.requires )
				if( strcmp( z.data.an{x}.requires{y}.name, oldName ) )
					z.data.an{x}.requires(y) = [];
					msg( sprintf( 'analysis unit ''%s'' required the unit you just deleted.\nthis requirement has been removed', z.data.an{x}.name ), 'msgbox' );
				end
			end
		end
	
		% first, assume all analysis units have callback functions...
		set( z.handle.widget.bottom.go, 'enable', 'on' );
		for x = 1 : length( z.data.an )
			if( ~length( z.data.an{x}.fun ) )
				% but if at least one unit doesn't have a CB function, disable execute button
				set( z.handle.widget.bottom.go, 'enable', 'off' );
			end	
		end


	end

	guidata( fig, z ) 

	updateUnitListDisplay( fig );

	% needed to resize analysis listbox height
	positionGUIElements( z );






% ===================================
function [] = newAnalysisCB( inarg1, inarg2, fig )

	z = guidata( fig );

	anNameCell = inputdlg( 'Enter name of new analysis unit:', '', 1 );
	if( length(anNameCell) ~= 1 )
		msg( 'new unit cancelled', 'status' );
		return; % user pressed cancel button
	end
	anName = cell2str( anNameCell );
	if( length( anName ) < 1 )
		msg( 'new unit cancelled', 'status' );
		return; % user pressed "ok" but didn't enter any text
	end
	msg( sprintf( 'analysis unit ''%s'' added', anName ), 'status' );

	% we have something to save now
	set( z.handle.widget.bottom.save, 'enable', 'on' );
	set( z.handle.widget.bottom.includeResults, 'enable', 'on' );

	% the new unit will be at the end of the list
	num = length( z.data.an ) + 1;

	if( analysisNameExists( anName, fig, [] ) == 1 ) 
		msg( 'the name you entered is already in use', 'msgbox' );
		return;
	else
		z.data.an{num} = newAnalysisUnit( anName );
	end

	% don't shift the focus of the currently selecetd unit, 
	% unless the new unit we just created is now the only one. 
	if( num == 1 ) 
		set( z.handle.widget.mid.jobSpecListbox, 'value', num );
	end

	% now we have something to delete
	set( z.handle.widget.top.deletean, 'enable', 'on' );
	
	% now we can apply the selected analysis funciton from right pane to the current an unit
	set( z.handle.widget.bottom.anAdd, 'enable', 'on' );

	% allow reordering of an unit list if there's at least 2
	if( length( z.data.an ) > 1 )
		set( z.handle.widget.bottom.moveup, 'enable', 'on' );
		set( z.handle.widget.bottom.movedown, 'enable', 'on' );
	end

	% since we have at least one analysis unit, we have one selected. 
	% now, if there's also some items selected in the left window, enable
	% the filt include/exclude buttons
	if( length( z.data.left.selected ) )
		set( z.handle.widget.bottom.addexclude, 'enable', 'on' );
		set( z.handle.widget.bottom.addinclude, 'enable', 'on' );
	else
		set( z.handle.widget.bottom.addexclude, 'enable', 'off' );
		set( z.handle.widget.bottom.addinclude, 'enable', 'off' );
	end

	% disable executwe button because we know that at least one analysis unit in
	% the list (the one we just added) doesn't have a callback function specified yet.
	set( z.handle.widget.bottom.go, 'enable', 'off' );


	% upload changes to guidata	
	guidata( fig, z ) 

	% needed to regenerate the analysis list options
	updateUnitListDisplay( fig );

	% need to resize analysis listbox height
	positionGUIElements( z );

	% bizarre matlab GUI bug- for some reason calling this again here fixes
	% a transient formatting problem which appears after adding the first unit
	resizeMainFigureCB( '', '', fig )




% ============================================
% generates human readable string describing the input filt. 
% returns a cell array of strings describing the input struct,
% which is intended to be a filt for this purpose, but this function will work with
% any struct. in fact ### we would want to pull this out into a static utility function
% if anybody else need it. 
% each string is one line of the filt such as: 
% filt.exclude.any.subject_id = {'tmp_*','01zin79271','08tgs78071'}
% filt.exclude.any.session_id = [1873 1984  1523:1576]
%
function [s] = getFilterList( f )

	s = {};
	if( length( f ) < 1 )
		return;
	end

	fn = fieldnames( f );
	howManyTotal = 0;

	for x = 1 : length( fn )
		fieldsIJustGot = {};

		nextField = getfield( f, fn{x} );
		if( ~isstruct(nextField) )
			prefChar = ' = ';
			if( iscell( nextField ) )
				temp = '{';
				for z = 1 : length( nextField )
					temp = sprintf( '%s''%s''', temp, nextField{x} );
					% add a comma following each element except the last in the cell array
					if( z < length(nextField) )
						temp = sprintf( '%s,', temp );
					end
				end
				temp = sprintf( '%s}', temp );
				fieldsIJustGot{1} = temp;
			else
				if( isnumeric( nextField ) )
					temp = '[';
					for z = 1 : length( nextField )
						temp = sprintf( '%s%d', temp, nextField(z) );
						% add a space following each element except the last in the cell array
						if( z < length(nextField) )
							temp = sprintf( '%s ', temp );
						end
					end
					temp = sprintf( '%s]', temp );
					fieldsIJustGot{1} = temp;
				else
					fieldsIJustGot{1} = '?DF34 data type not accounted for?';
				end
			end
	
		else
			prefChar = '.';
			% recursive function call
			fieldsIJustGot = getFilterList( nextField );
		end

		for y = 1 : length( fieldsIJustGot )
			howManyTotal = howManyTotal + 1; % initial value is 0
			s{howManyTotal} = sprintf( '%s%s%s', fn{x}, prefChar, fieldsIJustGot{y} );
		end

	end






% ============================================
% create a blank analysis unit
function [an] = newAnalysisUnit( name )

	% defaults
	an.type = '';
	an.fun = [];
	an.requires = {};
	an.name = name;
	an.params.filt = [];
	an.results = {};

	
	if( 0 ) % DEBUG 
		an.name = name;
		an.params.x1 = 123;
		an.params.y2 = 'blah blah';
		an.params.z3 = 12.2333;
		an.params.abc = [1 2.2 3];
		an.params.def = [1 2 3
							4 5 6
							7 8 9];
		an.params.col_vec = [3 5 23 7 4 2 1]';
		e.x1 = 342;
		e.x2 = '4534jlkj';
		e.x3 = [1 5 6];
		f.a1 = 23;
		f.a2 = e;
		an.params.sss = f;
		cell{1} = 'ghgfhf';
		cell{2} = 'gkfkh';
		cell{3} = '45353s';
		an.params.xyz = cell;
	end

	if( 0 ) %  DEBUG
		saveP = an.params;
		an = generateTestUnit;
		an.params = saveP;
		an.params.filt = [];
		an.results = {};
		an.type = '';
	end






% ===================================
function [] = addToFilterCB( inarg1, inarg2, fig )
	z = guidata( fig ); 

	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );

	if( strcmp( get( gcbo, 'tag' ), 'addinclude' ) )
		isInclude = 1;
	else
		isInclude = 0;
	end

	sel = z.data.left.selected;
	filtAddition = [];
	for x = 1 : size( sel, 1 )
		if( sel(x,1) > 0 )
			if( isInclude )
				filtAddition.include.all.experiment_id = [sel(x,1)];
			else
				filtAddition.exclude.all.experiment_id = [sel(x,1)];
			end
		end
		if( sel(x,2) > 0 )
			if( isInclude )
				filtAddition.include.all.form_id = [sel(x,2)];
			else
				filtAddition.exclude.all.form_id = [sel(x,2)];
			end
		end
		if( sel(x,3) > 0 )
			if( isInclude )
				filtAddition.include.all.question_id = [sel(x,3)];
			else
				filtAddition.exclude.all.question_id = [sel(x,3)];
			end
		end
		z.data.an{whichAnUnit}.params.filt = add2filt( z.data.an{whichAnUnit}.params.filt, filtAddition );
	end

	guidata( fig, z );

	% create a new ui tree because the existing one is now displaying stale filt info
	regenerateTreeForSelectedJob( '', '', fig );





% ====================================
% this string gets assigned to the listbox for the analysis units.
% each list display item is separated by a | character. 
function [listString] = getAnalysisUnitListString( fig )
	z = guidata( fig );
	listString = '';

	if( length( z.data.an ) < 1 )
		return; 
	end

	listString = ['1) ' z.data.an{1}.name]; 
	for x = 2 : length( z.data.an )
		listString = sprintf( '%s|%d) %s', listString, x, z.data.an{x}.name );
	end




% ====================================
% regenerates the list based on the available analysis units. 
function [] = updateUnitListDisplay(fig)

	z = guidata( fig );

	if( length( z.data.an ) < 1 )
		set( z.handle.widget.mid.jobSpecListbox, 'string', '' );
	
		% only change the current selection if the list is empty- must set to 0.
		set( z.handle.widget.mid.jobSpecListbox, 'value', 0 );
	else
		set( z.handle.widget.mid.jobSpecListbox, 'string', getAnalysisUnitListString( fig ) );
	end


	guidata(fig,z);







% ======================================
function [] = selectExperimentListCB( a1, a2, fig )

	z = guidata( fig ); 
	v = get(gcbo,'value'); % v is the list element number (not the experiment ID)
	expNum = z.data.meta.expList{1}(v); % here is the actual experiment ID
	
	% destroy any previous experiment selection uitree
	% isn't there an easier way to check for the existence of a
	% previous graphics object with this handle than using
	% isfield at every level? 
	h = z.handle;
	if( isfield( h, 'menu' ) )
		m = h.menu;
		if( isfield( m, 'leftTree' ) )
			delete( z.handle.menu.leftTree );
			z.handle.menu = rmfield( z.handle.menu, 'leftTree' );
			msg( 'deleting previous uitree', 'status' );
		end
	end

	msg( sprintf( 'loading', expNum, cell2str(z.data.meta.expList{2}(v)) ), 'status' );

	% nothing is selected any more
	z.data.left.selected = [];
	set( z.handle.widget.bottom.addexclude, 'enable', 'off' );
	set( z.handle.widget.bottom.addinclude, 'enable', 'off' );
		

	% create ui tree for the selected experiment
	z.handle.menu.leftTree = ens_getLeftTree( expNum, z.handle.window.main );

	% you can't set the parent of the tree directly, you have to go through its container.
	z.handle.container.leftTree = get(z.handle.menu.leftTree, 'uicontainer');
	set(z.handle.container.leftTree, 'parent', z.handle.panel.left );

	msg( sprintf( 'loaded experiment %d: %s', expNum, cell2str(z.data.meta.expList{2}(v)) ), 'status' );

	guidata( fig, z );




% ========================================================================
% general messaging fucntion supports types:
% status - writes to the status/info bar at bottom of GUI
% msgbox - popup message
% errbox - popup error message
% default - fprintf to matlab console output
function [] = msg( msg, outputType )

	if( ~exist( 'outputType', 'var' ) )
		outputType = 'default';
	end
		

	if( strcmp( outputType, 'errbox' ) )
		errordlg( msg, 'BAD' );
		return;
	end 

	if( strcmp( outputType, 'msgbox' ) )
		msgbox( msg, 'I would like you to know this:' );
		return;
	end
	
	if( strcmp( outputType, 'status' ) )
		fig = getOurGUIHandle;

		% if it didn't work, something is really wrong.
		if( ~length(fig) )
			return;
		end

		% check to see if this is too long for the status bar?
		% right now it just stops at the end. maybe coudl at least replace with '...'
		z = guidata( fig );
		set( z.handle.text.status, 'string', msg );
		return;
	end

	if( strcmp( outputType, 'homing_pigeon' ) )
		% ### need to register homing pigeon service with matlab
	end

	% default
	fprintf( '%s\n', msg );
	
	


% ======================================
% this is for functions that don't have direct access to our figure.
% this is much safer than using gcf()
function [h] = getOurGUIHandle()

	h = [];
	c = get(0,'children');
	for x = 1 : length(c)
		if( strcmp( get(c(x),'tag'), 'z.handle.window.main' ) )
			h = c(x);
			return;
		end
	end

	fprintf( 'ERROR: couldn''t find main GUI figure handle\n\n' );





% ==============================================
% return the first node in the uitree for the selected analysis unit
function [tree] = getMiddleTree( an, fig )
	
	z = guidata( fig );

	% uitreenode( value, description, icon, leaf )
	% nodes values are assigned as: <variable name>.<variable value>
	% we have to keep track of these name.value pairs separately
	% from what is displayed to the user, which might not necessarily
	% always want to be the same.
	value = sprintf( 'name.%s', an.name );
	rootNode = uitreenode( value, an.name, [], false ); 

	% attach the rootNode to the uitree and set the expansion callback.
	tree = uitree(	'Root', rootNode, ...
								'ExpandFcn', {@expandAnalysisUnitTreeCB, fig} );
	set( tree, 'NodeSelectedCallback', {@analysisFieldSelectionCB, fig} );
	set( tree, 'MultipleSelectionEnabled', 0 );
	% position and size are relative to the parent uicontainer or figure. 
	set(tree, 'Units', 'normalized', 'position', [0 0 1 1] );




% =======================================================
function analysisFieldSelectionCB( tree, value, fig )
	z = guidata( fig );

	nodes = tree.SelectedNodes;

	% ### multiple selection is disabled. 
	% things would get tricky and confusing fast if multiple analysis unit selection were allowed
	if( length( nodes ) ~= 1 )
		msg( sprintf( '%d nodes selected in analysis unit uitree (should be 1)', length(nodes)), 'errbox' );
		msg( 'aborted', 'status' );
		return;
	end

	% node value which is created by us comes in the form <name>.<value>
	% there are some nodes that have multiple values: <name>.<value1>.<value2>.<value3>...
	% but we are only splitting at the first "." right now. anything to the right
	% of the first . is considered the entire "value". later when each node is dealt
	% with individually based on its name, those special sections will parse out
	% the the value to each of its components if necessary. see: params for the main example  
	tag = get( nodes(1), 'value' );
	p = strfind( tag, '.' );
	if( p <= 0 )
		msg( sprintf( 'bad node value: %s\n', tag ), 'errbox' );
		return;
	end

	% we could pass these directly as params to updateFieldEditRegion but it's a good idea to store 
   z.data.center.vName = tag(1:p-1);
   z.data.center.vValue = tag( p+1:end ); % might contain multiple sub-values, each separated by a "."

	if( 0 ) % DEBUG
		msg( sprintf( 'selected field name -->%s<--  value-->%s<--', z.data.center.vName, z.data.center.vValue ), 'status' );
	end

	% save changes
	guidata( fig, z );


	% this will create all the necessary input fields in the center
	% bottom panel for editing analysis field values
	updateFieldEditRegion( fig );




% =========================================================
function [] = disableAllfieldEditWidgets( fig )
	z = guidata( fig );

	set( z.handle.text.fieldName, 'string', '' );

	if( strcmp( get(  z.handle.widget.mid.fieldEditPopupmenu, 'visible' ), 'on' ) )
		 set( z.handle.widget.mid.fieldEditPopupmenu, 'visible', 'off' );
		 set( z.handle.widget.mid.fieldEditPopupmenu, 'enable', 'off' );
	end
	if( strcmp( get(  z.handle.widget.mid.fieldEditString, 'visible' ), 'on' ) )
		 set( z.handle.widget.mid.fieldEditString, 'visible', 'off' );
		 set( z.handle.widget.mid.fieldEditString, 'enable', 'off' );
	end

	if( strcmp( get(  z.handle.widget.mid.filtDeleteButton, 'visible' ), 'on' ) )
		 set( z.handle.widget.mid.filtDeleteButton, 'visible', 'off' );
		 set( z.handle.widget.mid.filtDeleteButton, 'enable', 'off' );
	end






% =========================================================
% looks at the global vNAme and vValue to see which node is currently selected.
% based on each special case, this function determines which editing options
% should be displayed to the user at the center-bottom panel, if any.
% for example if vName is "name", we enable the generic text field
% to allow the user to change the name of the selected analysis unit. 
function [] = updateFieldEditRegion( fig )
	z = guidata( fig );

	% these store the current node selection
	vName = z.data.center.vName;
	vValue = z.data.center.vValue;

	% first make any previously used field invisible/disabled
	disableAllfieldEditWidgets( fig );

	% if any of these is selected, we want to leave the edit region disabled- 
	% there's nothing that the user can change about any of these fields
	if( 	strcmp( vName, 'fun' ) | ...
			strcmp( vName, 'results' ) |	strcmp( vName, 'filtItem' ) )
		return;
	end

	if( strcmp( vName, 'filt' ) )
		set( z.handle.widget.mid.filtDeleteButton, 'visible', 'on' ); 
		set( z.handle.widget.mid.filtDeleteButton, 'enable', 'on' ); 
		set( z.handle.text.fieldName, 'string', 'filt' );
		return;
	end
	

	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	if( whichAnUnit < 1 | whichAnUnit > length( z.data.an ) )
		msg( 'error VN57', 'errbox' );	
		return;
	end	


	% we need to use the popup list for this
	if( strcmp( vName, 'requires' ) )

		listString = getAnalysisUnitListString( fig );
		set( z.handle.widget.mid.fieldEditPopupmenu, 'string', listString );
		set( z.handle.widget.mid.fieldEditPopupmenu, 'visible', 'on' ); 
		set( z.handle.widget.mid.fieldEditPopupmenu, 'enable', 'on' ); 
	
		% go through this whole thing just to find out which 
		% unit the selected unit currently requires, if any. 
		% if so, we want to set the initial value of the popup window
		% to the required unit for clarity and smooth operation.
		% sometimes the popup window is glitchy as it is. 
		whichRequired = 1;
		if( length( z.data.an{whichAnUnit}.requires ) ) 
			for t = 1 : length( z.data.an )
				if( strcmp( z.data.an{whichAnUnit}.requires{1}.name, z.data.an{t}.name ) )
					whichRequired = t;
					break;
				end
			end
		end

		set( z.handle.widget.mid.fieldEditPopupmenu, 'value', whichRequired );
		set( z.handle.text.fieldName, 'string', vName );

		return;
	end

	if( strcmp( vName, 'params' ) )

		if( length( vValue ) < 1 )
			% jg1
			% enable edit box to manually add fields to params struct.
			% commands are interpreted as if typed in the command line
			% following "an{this}.params." --- then the user might type for example
			% "t = 13;", restuing in "an{this}.params.t=13;"
			set( z.handle.widget.mid.fieldEditString, 'visible', 'on' ); 
			set( z.handle.widget.mid.fieldEditString, 'enable', 'on' ); 
			set( z.handle.widget.mid.fieldEditString, 'string', '' );
			set( z.handle.text.fieldName, 'string', 'cmd: params.' );
			return;	
		else
			% vValue might be a long list of nested stuct names like "s1.s2.s3.s4.selected_field_name'
			evalStr = sprintf( 'thisField = z.data.an{whichAnUnit}.params.%s;', vValue );
			eval( evalStr );
		end

		% nothing to modify here
		if( isstruct( thisField ) )
			return;
		end


		% we need to know the final field name for display purposes below, but we
		% don't have just the fnial name isolated by itself yet.
		% so we use this function as a convenience to parse the names for us, and
		% the final field name will be the last cell returned by the parse function
		names = parseNodeName( vValue );
		fieldName = names{ length( names ) }; 


		% most cases will use this
		set( z.handle.widget.mid.fieldEditString, 'visible', 'on' ); 
		set( z.handle.widget.mid.fieldEditString, 'enable', 'on' ); 


		[dataType, numDim] = getDataType( thisField );
		switch( dataType )

		% STRING
		case 's'
			set( z.handle.widget.mid.fieldEditString, 'string', thisField );


		% INT/FLOAT
		case {'d', 'f'}
			switch numDim	
			case 0 
				set( z.handle.widget.mid.fieldEditString, 'string', num2str( thisField ) );
			case 1 
				if( size(thisField,1) > 1 )
					set( z.handle.widget.mid.fieldEditString, 'string', ['[' num2str(thisField') ']'] );
					msg( 'NOTE: column vector displayed as row vector. if you save changes, it will become a row vector', 'status' );
				else
					set( z.handle.widget.mid.fieldEditString, 'string', ['[' num2str(thisField) ']'] );
				end
			otherwise	
				set( z.handle.widget.mid.fieldEditString, 'string', sprintf( '<%dd matrix>', numDim ) );
			end


		% CELL ARRAY 
		case {'c'}
			disableAllfieldEditWidgets( fig );	% this is kind of inefficient since we just enabled fieldEditString above
															% but there's no way around some repeated or redundant code with all these
			return;


		otherwise
			msg( 'error: ASL45', 'errbox' );
			return;
		end

		% if we made it here, then we're displaying something for the user to edit.
		% so change the text on the lower-left-center to the field name instead of "<no edit>"
		names = parseNodeName( vValue );
		set( z.handle.text.fieldName, 'string', names{ length(names) } );
		return;
	end


	% anything that makes it down here uses the simple text field 
	set( z.handle.widget.mid.fieldEditString, 'visible', 'on' ); 
	set( z.handle.widget.mid.fieldEditString, 'enable', 'on' ); 
	set( z.handle.widget.mid.fieldEditString, 'string', vValue );
	set( z.handle.text.fieldName, 'string', vName );






% ======================================================
% deletes the entire filter at once. 
% ### line-item deletion is not suported in v1.0
function [] = filtDeleteCB( var1, var2, fig )
	z = guidata( fig );

	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	newValue = get( gcbo(), 'value' );
	if( whichAnUnit < 1 | whichAnUnit > length( z.data.an ) )
		msg( 'eror FM101.3', 'errbox' );	
		return;
	end	

	z.data.an{whichAnUnit}.params.filt = [];
	msg( sprintf( 'filt for analysis unit %d has been deleted. see for yourself', whichAnUnit ), 'status' );

	guidata( fig, z );
	regenerateTreeForSelectedJob( '', '', fig );




% ======================================================
function [] = fieldEditPopupCB( var1, var2, fig )
	z = guidata( fig );

	% right now this only applies to selecting a new "requires" value-
	% need to deal with multiple possibilities for vName here
	% if other fields have popup menus. 

	% the current selection will always apply to one of the
	% analysis units in the list. it should never be possible to
	% get an out of bounds value here
	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	newValue = get( gcbo(), 'value' );
	if( whichAnUnit < 1 | whichAnUnit > length( z.data.an ) | newValue < 1 | newValue > length( z.data.an ) )
		msg( 'eror PY45', 'errbox' );	
		return;
	end	

	msg( sprintf( 'unit %s now requires unit %s', z.data.an{whichAnUnit}.name, z.data.an{newValue}.name ), 'status' ); 

	if( whichAnUnit == newValue )
		msg( 'warning: you have specified this unit''s output as its own input', 'msgbox' );
	end

	requiresStruct.name = z.data.an{newValue}.name;
	z.data.an{whichAnUnit}.requires{1} = requiresStruct;

	guidata( fig, z );
	regenerateTreeForSelectedJob( '', '', fig );




% ======================================================
function [res] = analysisNameExists( findThis, fig, exceptThese ) 
	z = guidata( fig );
	res = 0;

	if( length( z.data.an ) < 1 )
		return;
	end

	for x = 1 : length( z.data.an )
		if( length( find( exceptThese == x ) ) > 0 )
			continue; 
		end
	
		if( strcmp( z.data.an{x}.name, findThis ) )
			res = 1;
			return;
		end
	end







% ======================================================
function [] = fieldEditTextCB( var1, var2, fig )
	z = guidata( fig );

	
	vName = z.data.center.vName;
	vValue = z.data.center.vValue;

	newValue = get( gcbo(), 'string' );


	% the field chagne will always apply to one of the
	% analysis units in the list, so make sure we're in bounds
	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );
	if( whichAnUnit < 1 | whichAnUnit > length( z.data.an ) )
		msg( 'error TR32', 'errbox' );	
		return;
	end	


	if( strcmp( vName, 'name' ) )
		if( length( newValue ) < 1 )
			msg( 'can''t change name to blank', 'msgbox' );
			return;
		end
		% make sure name is not in use
		if( analysisNameExists( newValue, fig, whichAnUnit ) )
			msg( 'the name you entered is already in use', 'msgbox' );
			return;
		end
		z.data.an{whichAnUnit}.name = newValue;

		% if any existing unit required the unit whose name was just changed, 
		% change their requires field to match the new name to match the new name
		for x = 1 : length( z.data.an )
			if( ~length( z.data.an{x}.requires ) )
				continue;
			end
			if( strcmp( z.data.an{x}.requires{1}.name, vValue ) )
				msg( sprintf( 'unit %d (%s) requires the unit whose name you just changed.\nupdating %s''s requires field with the new name: %s', x, z.data.an{x}.name, z.data.an{x}.name, newValue ), 'msgbox' );
				z.data.an{x}.requires{1}.name = newValue;
			end
		end

	end


	if( strcmp( vName, 'type' ) )
		if( length( newValue ) < 1 )
			msg( 'can''t change type to blank', 'msgbox' );
			return;
		end
		z.data.an{whichAnUnit}.type = newValue;
	end
	

	if( strcmp( vName, 'params' ) )

		% jg2
		if( length( vValue ) <= 1 )
			% if here, then they are at the params level, and whatever they submitted
			% should be executed literally as an addition to params.
			%jg1 = 5
			if( length(newValue) <= 1 )
				msg( 'nothing entered', 'status' )
			else
				cmd = sprintf( 'z.data.an{%d}.params.%s;', whichAnUnit, newValue );
				try
					eval( cmd )
					msg( sprintf( 'executed: %s', cmd ), 'status' );
				catch 
					msg( sprintf( 'error executing: "%s"', cmd ), 'status' );
				end
			end
					
		else
			names = parseNodeName( vValue );

			finalFieldName = names{ length(names) };
			baseStructPath = 'params';
			for x = 1 : length( names )-1
				baseStructPath = sprintf( '%s.%s', baseStructPath, names{x} );
			end
			msg( sprintf( 'set param field %s.%s to: %s', baseStructPath, finalFieldName, newValue ), 'status' );
	
	
			% the newValue for this param comes to us as a string through the editable
			% text field, so we have to cast it to its true type.
			% the only way to find out its true type is to look at the previous
			% value of that param field
			evalStr = sprintf( 'oldFieldValue = getfield( z.data.an{whichAnUnit}.%s, ''%s'' );', baseStructPath, finalFieldName );
			eval( evalStr );
			[dataType, numDim] = getDataType( oldFieldValue );
			switch( dataType )
	
			case 's'	
				% no conversion necessary for newValue - leave as string
	
			case {'d','f'}
				% evalin is just like eval but lets you use the context of the workspace,
				% otherwise workspace variables aren't visible from here using eval()
				try
					newValue = evalin( 'base', newValue );
				catch
					msg( 'the string you entered didn''t evaluate to a valid number/vector/matrix', 'msgbox' );
					return;
				end
	
			otherwise
				msg( 'unknown data type err: MF321', 'errbox' );
			end
	
			% now assign newValue to the right place
			evalStr = sprintf( 'z.data.an{whichAnUnit}.%s.%s = newValue;', baseStructPath, finalFieldName );
			eval( evalStr );
				
		end % else - length(vValue) is <= 1
	end

	% save changes
	guidata( fig, z );

	% gnerates a new guitree
	regenerateTreeForSelectedJob( '', '', fig );




% separates the input string "<v1>.<v2>.<v3>.<v4>..."
% into a cell array: { 'v1', 'v2', 'v3', 'v4', ... }
function [v] = parseNodeName( s )

	p = strfind( s, '.' );
	if( length( p ) <= 0 )
		v{1} = s;
		return;
	end

	v{1} = s(1:p(1)-1);
	rem = s(p(1)+1:end);

	count = 2;
	while( 1 )
		p = strfind( rem, '.' );

		if( length(p) <= 0 )
			if( length( rem ) > 0 )
				v{count} = rem;
			end
			return; 
		else
			v{count} = rem( 1:p(1)-1 );
			rem = rem( p(1)+1:end );
			count = count + 1;
		end
	end
		



% ======================================================
% README:
% node names are kind of complicated because there are so many special
% cases- there's no single system or solution i could find which didn't
% have exceptions so i picked the system that made the most sense to me. 
% the node "value" is the part the the code uses to identify and modify
% values. generally the value is of the form <node name>.<node value>
% so for example the value of the node for "name" (analysis unit name)
% might be  value = "name.my_unit_1"
% however, things get more tricky when it comes to params, which 
% supports arbitrary nested structs. we use the node value here
% to store the path through the nested structs to the current node. 
% they all begin with "params", and one possible example might be:
% value = "params.abc.def.object34.object676"
% note that in this case the final element in that concatenated list of fields
% does not store the value of the field, as was the case with the 
% top-level node values (like name, function, requires, etc.)
% the top level root node for the params struct has value = "params."
%
% each of the special node types needs to be dealt with separately
% in 3 places: 1 is here, where the nodes are created when the user expands a
% section, and 2 when the user selects a node, we update the global vName and vValue,
% and 3 when the user saves changes to a field, we look at the stored vName and vValue  
% and deal with each different case accordingly.
function [nodes] = expandAnalysisUnitTreeCB( tree, value, fig )

	p = strfind( value, '.' );
	if( length( p ) < 1 )
		msg( sprintf( 'expandAnalysisUnitTreeCB: ERROR: bad node value: %s\n', value ), 'errbox' );
		nodes = [];
		return;
	end

	vName = value(1:p(1)-1);
	vValue = value(p(1)+1:end);

	if( 0 ) % DEBUG
		msg( sprintf( 'expanding analysis unit   %s.%s', vName, vValue ), 'status' );
	end

	z = guidata( fig );
	an = z.data.an{get(z.handle.widget.mid.jobSpecListbox,'value')}; 
	

	% if true then we're at the top level expansion
	if( strcmp( vName, 'name' ) )

		% this should never be possible
		if( ~strcmp( vValue, an.name ) )
			msg( sprintf( 'error X11: analysis unit names don''t match:  %s | %s', vValue, an.name ), 'errbox' );
			nodes = [];
			return;
		end

		val = sprintf( 'type.%s', an.type ); 
		if( length( an.type ) )
			desc = sprintf( 'type: %s', an.type );
		else
			desc = 'type: <empty> ';
		end
		nodes(1) = uitreenode( val, desc, [], true );
	
		if( length( an.fun ) )
			val = sprintf( 'fun.%s', func2str( an.fun ) ); 
			desc = sprintf( 'function: @%s', func2str( an.fun ) );
		else
			val = 'fun.';
			desc = 'function: <none specified>';
		end
		nodes(2) = uitreenode( val, desc, [], true );
	
		val = 'requires.';
		if( length( an.requires ) )
			% TODO: multiple cells in the reuiqres field aren't really supported in the gui right now. 
			% if an outside analysis unit is somehow ipmorted with multiple reuqires they will be displayed
			% here, but the user only has the ability to modify the first requires{1} cell
			desc = sprintf( 'requires: %s', an.requires{1}.name );
			for r = 2 : length( an.requires )
				desc = sprintf( '%s, %s', desc, an.requires{r}.name );
			end
		else
			desc = 'requires: <no input selected>';
		end
		nodes(3) = uitreenode( val, desc, [], true );
	
		val = sprintf( 'params.' ); 
		desc = sprintf( 'params' );
		nodes(4) = uitreenode( val, desc, [], false );

		val = sprintf( 'results.' ); 
		desc = sprintf( 'results: %d cells', length( an.results ) );
		nodes(5) = uitreenode( val, desc, [], true );

		return;
	end

	if( strcmp( vName, 'params' ) )

		names = parseNodeName( value );
		structPath = 'params';
		for x = 2 : length( names )
			structPath = sprintf( '%s.%s', structPath, names{x} );
		end
		evalStr = sprintf( 's = an.%s;', structPath );
		eval( evalStr );

		fn = fieldnames( s );
		for x = 1 : length(fn)
			if( strcmp( fn{x}, 'filt' ) )
				val = sprintf( 'filt.' ); 
				desc = sprintf( 'filt' );
				nodes(x) = uitreenode( val, desc, [], false );
			else 
				val = sprintf( '%s.%s', structPath, fn{x} ); 
				fieldValue = getfield( s, fn{x} );
				[dataType, numDimensions] = getDataType( fieldValue );
				isLeaf = true;
				switch dataType

				case 's'	
					desc = sprintf( '%s: ''%s''', fn{x}, fieldValue ); 

				case {'d','f'}
					if( numDimensions < 2 )
						if( numDimensions == 1 )
							if( size( fieldValue, 1 ) > 1 )
								desc = sprintf( '%s: [%s] (column)', fn{x}, num2str(fieldValue') ); 
							else
								desc = sprintf( '%s: [%s]', fn{x}, num2str(fieldValue) ); 
							end
						else
							desc = sprintf( '%s: %s', fn{x}, num2str(fieldValue) ); 
						end
					else
						desc = sprintf( '%s: %dd matrix', fn{x}, numDimensions ); 
					end


				case 'o'
					desc = sprintf( '%s', fn{x} ); 
					isLeaf = false;

				case 'c'
					desc = sprintf( '%s: %d cells', fn{x}, length(fieldValue) ); 

				otherwise
					desc = 'unknown field type';
					msg( 'unknown field type  BB29', 'errbox' );
				end

				nodes(x) = uitreenode( val, desc, [], isLeaf );
			end
		end

		return;
	end

	if( strcmp( vName, 'filt' ) )
	
		filtList = getFilterList( an.params.filt );
		if( length( filtList ) )
			for x = 1 : length( filtList )
				val = sprintf( 'filtItem.%d', x );
				desc = sprintf( '%s\n', filtList{x} );
				nodes(x) = uitreenode( val, desc, [], true );
			end
		else
			nodes = [];
		end

		return;
	end

	msg( sprintf( 'expand analysis unit: unknown node name: %s', vName ), 'errbox' );
	nodes = [];





%=================================================
% numDimensions:
% 0 is a single value
% 1 is a 1d vector
% 2 is a 2d matrix 
% N is an Nd matrix
%
% t:
% o = struct
% c = cell array
% d = int
% f = float
% s = string
% ? = ?
function [t, numDimensions] = getDataType( in )
	numDimensions = 0;

	if( isstruct( in ) )
		t = 'o';
		return;
	end
	
	if( iscell( in ) )
		t = 'c';
		return;
	end

	if( ischar( in ) )
		t = 's';
	else
		if( isnumeric( in ) )
			if( numel( in ) > 1 )
				numDimensions = length(size(in)); % same as ndims(in)
			end
			if( numDimensions == 2 )
				% at this point we know we have either a vector or a 2D matrix.
				% but the caller needs to know which one it is, so here we
				% decrement numDimensions from 2 to 1 for vectors. 
				if( size(in,1) ==1 | size(in,2) == 1 )
					numDimensions = numDimensions - 1;	
				end
			end
				
			if( ceil( in ) == in )
				t = 'd';
			else
				t = 'f';
			end
		else
			t = '?';
		end
	end



function [tree] = ens_getLeftTree( experimentNumbers, figHandle )
	
	% if you are implementing support for multiple experiments, 
	% you must search all code for: SME
	% this means something about that code needs to be
	% changed in order to enable multiple experiment support.
	% these markers aren't necessarily the only changes that need to be made
	if( length( experimentNumbers ) ~= 1 )
		% ### SME
		fprintf( 'filterSelectionUI: ERROR: multiple experiments not supported\n' );
		tree = [];
		return;
	end

	z = guidata( figHandle );


	expMetaData = mysql_extract_metadata( ...
											'table', 'experiment', ...
											'conn_id', z.config.conn_id, ...
											'experiment_id', experimentNumbers ); 

	% value, description, icon, leaf
	% nodes values are assigned as: level.row
	% so this root node is at experiment level and it's the 1 and only. 
	rootNode = uitreenode( 'experiment.1', expMetaData.experiment_title, [], false ); 
	set( rootNode, 'userdata', experimentNumbers); % ## SME

	% attach the rootNode to the uitree and set the expansion callback.
	tree = uitree(	'Root', rootNode, ...
								'ExpandFcn', {@expandGUITreeCB, expMetaData} );
	set( tree, 'NodeSelectedCallback', {@leftTreeSelectionCB, figHandle} );
	set( tree, 'MultipleSelectionEnabled', true );
	% position and size are relative to the parent uicontainer or figure. 
	set(tree, 'Units', 'normalized', 'position', [0 0 1 1] );




% =======================================================
function leftTreeSelectionCB( tree, value, figHandle )

	nodes = tree.SelectedNodes;

	if( 0 ) % DEBUG
		fprintf( '%d items selected\n', length(nodes) );
	end

	% data for selections keeps track of experiment number,
	% form number, and question number. 
	% a selection can happen at any level, so a selection 
	% line item (one node) may be an experiment, a form,
	% or a particular question. 
	% the selectedData matrix stores all mixed types of line
	% items, one per row, with cols 1,2,3 storing experiment,form,question
	% ids in that order. selections at the form level will leave
	% the 3rd column blank. example of a selectedData matrix: 
	%	[12		45		123
	%	12			0		0
	%	12			45		124
	%	12			46		0]
	selectedData = zeros(length(nodes),3);

	for x = 1 : length( nodes )

		% parse the value field of this node to find out
		% if this is an experiment, form, or question
		v = get( nodes(x), 'value' );
		p = strfind( v, '.' );
		type = v(1:p-1);

		% we could add error checking for p < 0 but this should never happen

		if( strcmp( type, 'experiment' ) )
			selectedData(x,:) = [get(nodes(x),'userdata') 0 0];
		end

		if( strcmp( type, 'form' ) )
			% find out experiment id by looking at parent
			e = get(get(nodes(x),'parent'),'userdata');
			selectedData(x,:) = [e get(nodes(x),'userdata') 0];
		end
		
		if( strcmp( type, 'question' ) )
			% find out experiment id by looking at parent
			p1 = get(nodes(x),'parent');
			p2 = get(p1,'parent');
			f = get(p1,'userdata');	
			e = get(p2,'userdata');	
			q = get(nodes(x),'userdata');
			selectedData(x,:) = [e f q];
		end

		if( 0 ) % DEBUG
			fprintf( 'selected %d: [%d %d %d]\n', x, selectedData(x,1),selectedData(x,2),selectedData(x,3) );
			fprintf( '%d name: %s   value: %s  userdata: %d\n', ...
							x, get(nodes(x),'name'), get(nodes(x), 'value'), get(nodes(x),'userdata') );
		end
	end

	% update global (figure) app data struct
	z = guidata( figHandle );

	% toggle enable submit button for the left panel
	% i don't think there's any way to unselect all
	% and even if there were, i don' tthink it would 
	% call this selection CB, but handle the 0 case just in case
	
	% there has to be an existing analysis unit selected also
	whichAnUnit = get( z.handle.widget.mid.jobSpecListbox, 'value' );

	
	if( length(nodes) > 0 & whichAnUnit > 0 )
		set( z.handle.widget.bottom.addexclude, 'enable', 'on' );
		set( z.handle.widget.bottom.addinclude, 'enable', 'on' );
	else
		set( z.handle.widget.bottom.addexclude, 'enable', 'off' );
		set( z.handle.widget.bottom.addinclude, 'enable', 'off' );
	end

	z.data.left.selected = selectedData;
	guidata( figHandle, z );






% ======================================================
function [nodes] = expandGUITreeCB( tree, value, expMetaData )


	p = strfind( value, '.' );
	if( p <= 0 )
		fprintf( 'expandGUITreeCB: ERROR: bad node value: %s\n', value );
		nodes = [];
		return;
	end

	level = value(1:p-1);
	row = str2num( value( p+1:end) );

	if( 0 ) % DEBUG
		fprintf( 'expanding level %s, row %d\n', level, row );
	end
	
	if( strcmp( level, 'experiment' ) )
		numChild = length( expMetaData.form );
		nextLevel = 'form';
	end
	if( strcmp( level, 'form' ) )
		numChild = length( expMetaData.form(row).question );
		nextLevel = 'question';
	end
	
	% create a uitree node for each of the items in the next expansion
	for x = 1 : numChild	


		if( strcmp( level, 'experiment' ) )
			desc = sprintf( '%s (%d)', ...
 				expMetaData.form(x).form_name, ...
 				expMetaData.form(x).form_id );

			% if it has questions associated with it, then it's not a leaf
			if( length( expMetaData.form(x).question ) > 0 )
				isLeaf = false;
			else
				isLeaf = true;
			end
			record_id = expMetaData.form(x).form_id;
		end


		if( strcmp( level, 'form' ) )
			desc = sprintf( '%s (%s)', ...
 				expMetaData.form(row).question(x).question_text, ...
 				sprintf( '%.04f',expMetaData.form(row).question(x).compqid ) );
			isLeaf = true;
			record_id = expMetaData.form(row).question(x).question_id;
		end
	
		% matlab automatically attaches this list of nodes to the
		% parent in the uitree. 
		nodes(x) = uitreenode( sprintf( '%s.%d',nextLevel,x ), ...
						desc, [], isLeaf );
		set(nodes(x),'userdata',record_id);

	end
	






















