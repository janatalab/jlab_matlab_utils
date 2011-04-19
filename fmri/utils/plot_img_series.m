function plot_img_series(imgInfo,params)
% plot_img_series(imgInfo, params)
%
% Uses slice_overlay to plot multiple pages of slice data, with the
% possibility of multiple data volumes per page.
%
% imgInfo is equivalent to the SO structure that is expected by
% slice_overlay.
%
% Plotting parameters
% params.maxImgPerPage - maximum number of images per page [default=4]
%   .start_fignum - starting figure number [1]
%   .numCol - number of columns [2]
%   .rowTitleAreaHeight - proportion of row height devoted to plot titles [0.025]
%   .SEPARATE_FIGS - place each page in its own figure window [0]
%   .PRINT_FIGS - print the figure [0]
%   .figfname - path to which to print the figure ['']

% December 2010, Petr Janata

slice_overlay_defaults;

% See if we have an image transformation structure specified
if isfield(params,'tp')
	tp = params.tp;
end

% Check to see if we have necessary variables
try maxImgPerPage = params.maxImgPerPage; catch maxImgPerPage = 4; end
try SEPARATE_FIGS = params.SEPARATE_FIGS; catch SEPARATE_FIGS = 0; end
try start_fignum = params.start_fignum; catch start_fignum = 1; end
try numCol = params.numCol; catch numCol = 2; end
try rowTitleAreaHeight = params.rowTitleAreaHeight; catch rowTitleAreaHeight = 0.025; end
try colScaleFactor = params.colScaleFactor; catch colScaleFactor = 0.95; end
try PRINT_FIGS = params.PRINT_FIGS; catch PRINT_FIGS = 0; end
try figfname = params.figfname; catch figfname = ''; end

numRow = ceil(maxImgPerPage/numCol);

% Initialize plot parameters
clear global SO
global SO

% Figure out how many ICs we need to plot
numPlots = length(imgInfo);

numPages = ceil(numPlots/maxImgPerPage);
	
fignum = start_fignum;
for ipage = 1:numPages
	
	% Calculate offset into imgInfo
	pageOffset = (ipage-1)*maxImgPerPage+1;
	
	% Figure out indices in imgInfo that we'll be using
	pageImgIdxs = pageOffset:min(pageOffset+maxImgPerPage-1,numPlots);
		
	% Create a new figure and jump through some hoops to size it
	% correctly so that slice_overlay will work
	if ipage > 1 && SEPARATE_FIGS
		fignum = fignum+1;
	end
		
	% Position the figure appropriately
	figure(fignum), clf
	set(gcf,'PaperPosition', [0 0 1 1], 'PaperUnits','normalized'); % create a full page
	set(gcf,'PaperUnits', 'points', 'Units','points');  % convert to measure we can use to set the position with
	set(gcf,'Position', [0 0 get(gcf,'PaperSize')]); % move the figure window on the screen and set it to correct size
	set(gcf,'Units','pixels'); % put things in pixels which is what slice_overlay wants
		
	SO.figure = figure(fignum);
	set(SO.figure,'DefaultLineLineWidth', DEFAULT_CONTOUR_WIDTH)
	SO.labels.format = '%+d';
	SO.labels.size = 0.1;
	SO.transform = tp.transform_types{tp.transform_idx};
		
	for iplot = 1:length(pageImgIdxs)
		
		% Copy the img stack
		SO.img = imgInfo(pageImgIdxs(iplot)).img;
		
		try SO.cbar = imgInfo(pageImgIdxs(iplot)).cbar; catch SO.cbar = []; end
		try SO.contours = imgInfo(pageImgIdxs(iplot)).contours; catch SO.contours = []; end
		try SO.ticks = imgInfo(pageImgIdxs(iplot)).ticks; catch SO.ticks.show = 0; end		
    try SO.white_background = imgInfo(pageImgIdxs(iplot)).white_background; catch SO.white_background = 0; end
   
    if SO.white_background
      SO.labels.colour = [0 0 0];
    end
    
		% Specify the plot area
		rowidx = ceil(iplot/numCol);
		colidx = mod(iplot-1,numRow)+1;
			
		row_height = 1/numRow;
		height = row_height-rowTitleAreaHeight;
			
		col_margin = ((1-colScaleFactor)/numCol)/2;
		col_width = 1/numCol;
		width = col_width-col_margin*2;
			
		bottom = 1-rowidx*row_height;
		left = (colidx-1)/numCol+col_margin;
			
		% Make the plot
		if ~isempty(tp.non_contig_slices)
			SO.slices = tp.non_contig_slices;
		else
			SO.slices = ...
				tp.slice_ranges(tp.transform_idx,1):sign(diff(tp.slice_ranges(tp.transform_idx,:)))*tp.slice_skip:tp.slice_ranges(tp.transform_idx,2);
		end
		slice_overlay('checkso')
			
		SO.refreshf = 1;
		SO.clf = 0;
		SO.area.position = [left bottom width height];
		SO.area.units = 'normalized';
		slice_overlay('display')  % do the plot
			
		% Now we need to add a title
		try 
			title_str = imgInfo(pageImgIdxs(iplot)).title.text;
		catch
			title_str = '';
		end
		
		try 
			format_args = imgInfo(pageImgIdxs(iplot)).title.format_args;
		catch 
			format_args = {};
		end
			
		% Create an axes in the title area
		tah = axes('position', [left bottom+height width rowTitleAreaHeight], ...
			'visible','off');
			
		if ~isempty(title_str)
			t = text(0.5, 0.25, title_str, format_args{:});
		end
	end % for iplot=
		
	% Print the figure out to a file
	if PRINT_FIGS && ~isempty(figfname)
		if ipage == 1
			fprintf('Printing plots to %s\n', figfname);
			append_str = '';
		else
			append_str = '-append';
		end
		set(gcf,'PaperPosition',[0 0 get(gcf,'PaperSize')]);
		print(figfname, '-dpsc', append_str);
	end
end % for ipage
