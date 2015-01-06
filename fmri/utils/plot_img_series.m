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
%   .PLOT_INDIV_SLICES - plots each slice to its own jpeg file in a
%   directory 

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
try PLOT_INDIV_SLICES = params.PLOT_INDIV_SLICES; catch PLOT_INDIV_SLICES = 0; end
try PAPER_TYPE = params.paperType; catch PAPER_TYPE = 'usletter'; end

numRow = ceil(maxImgPerPage/numCol);

global_tp = tp;

% Initialize plot parameters
clear global SO
global SO

if maxImgPerPage == 1
  numCol = 1;
end

% Figure out how many total plots we have
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
  set(gcf,'PaperType', PAPER_TYPE);
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
    currIdx = pageImgIdxs(iplot);
		
		% Copy the img stack
		SO.img = imgInfo(currIdx).img;
		
		try SO.cbar = imgInfo(currIdx).cbar; catch SO.cbar = []; end
		try SO.contours = imgInfo(currIdx).contours; catch SO.contours = []; end
		try SO.ticks = imgInfo(currIdx).ticks; catch SO.ticks.show = 0; end		
    try SO.white_background = imgInfo(currIdx).white_background; catch SO.white_background = 0; end
   
    if SO.white_background
      SO.labels.colour = [0 0 0];
    end
    
		% Specify the plot area
		rowidx = ceil(iplot/numCol);
		colidx = mod(iplot-1,numCol)+1;
			
		row_height = 1/numRow;
		height = row_height-rowTitleAreaHeight;
			
		col_margin = ((1-colScaleFactor)/numCol)/2;
		col_width = 1/numCol;
		width = col_width-col_margin*2;
			
		bottom = 1-rowidx*row_height;
		left = (colidx-1)/numCol+col_margin;
			
		% Make the plot
    if isfield(imgInfo(currIdx), 'tp')
      tp = imgInfo(currIdx).tp;
    else
      tp = global_tp;
    end
    
    if ~isempty(tp.non_contig_slices)
      SO.slices = tp.non_contig_slices;
    else
      if length(tp.slice_ranges) > 1
        SO.slices = ...
          tp.slice_ranges(tp.transform_idx,1):sign(diff(tp.slice_ranges(tp.transform_idx,:)))*tp.slice_skip:tp.slice_ranges(tp.transform_idx,2);
      else
        SO.slices = tp.slice_ranges;
      end
    end
    
		slice_overlay('checkso')
			
		SO.refreshf = 1;
		SO.clf = 0;
		SO.area.position = [left bottom width height];
		SO.area.units = 'normalized';
		slice_overlay('display')  % do the plot
			
		% Now we need to add a title
		try 
			title_str = imgInfo(currIdx).title.text;
		catch
			title_str = '';
		end
		
		try 
			format_args = imgInfo(currIdx).title.format_args;
		catch 
			format_args = {};
		end
			
		if ~isempty(title_str)
      % Create an axes in the title area
      tah = axes('position', [left bottom+height width rowTitleAreaHeight], ...
			'visible','off');
			
			t = text(0.5, 0.25, title_str, format_args{:});
		end
	end % for iplot=
		
  % If we are plotting individual slices and have only one child axes, then
  % try to crop the figure to the axes
  children = get(gcf,'children');
  if length(children) == 1
    % Get the position of the child and simply try to set the paper
    % position to that
%    childPos = get(children,'Position');
%    set(gcf,'Position',childPos)
%    set(gcf,'PaperPositionMode','auto')
%    set(children,'Position',[0 0 1 1], 'units', 'normalized')
    
    % To make sure that tick labels show up, set their color to black
    set(children,'xcolor', [0 0 0], 'ycolor', [0 0 0])
  end
  
	% Print the figure out to a file
  if PLOT_INDIV_SLICES && exist(imgInfo(currIdx).fpath,'dir')
    transform_str = tp.transform_types{tp.transform_idx};

    [fpath,fstub] = fileparts(imgInfo(currIdx).fpath);
    switch transform_str
      case 'coronal'
        slice_plane = 'y';
      case 'axial'
        slice_plane = 'z';
      otherwise
        slice_plane = 'x';
    end
    figfname = fullfile(imgInfo(currIdx).fpath, sprintf('%s_%s=%dmm.png',fstub, slice_plane,tp.slice_ranges));
    append_str = '';
    ftype_str = '-dpng';
  end 
  
	if PRINT_FIGS && ~isempty(figfname)
    if ~PLOT_INDIV_SLICES
      ftype_str = '-dpsc';
      if ipage == 1
        fprintf('Printing plots to %s\n', figfname);
        append_str = '';
      else
        append_str = '-append';
      end
    end
		set(gcf,'PaperPosition',[0 0 get(gcf,'PaperSize')]);
		print(figfname, ftype_str, append_str);
	end
end % for ipage
