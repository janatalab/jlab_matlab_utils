function add_errorbars(h,errordata,linecolor,orientation)
%
% add_errorbars(h,errordata);
%
% Adds error bars to bar plots.  The bar information is in h which is a vector
% of object handles returned by bar().  Error data is the error data to be
% plotted.

% 02/09/05 PJ -- Added a bit more dimension error checking, and object type
% checking
% 11/03/06 PJ -- Added handling of Matlab 7 barseries (hggroup) objects
% 11/11/06 PJ -- Fixed handling of single barseries object
% 05/06/07 PJ -- Added support for horizontal bargraphs

ncond = length(h);

if nargin < 3
  linecolor = '';
end

if nargin < 4
  orientation = 'vertical';
end

obj_type =  get(h(1),'Type');
switch obj_type
  case 'patch'
    
  case 'hggroup'
    % The patch object is actually buried deeper down. It is a child of the
    % hggroup object
    
    new_kids = get(h,'Children');
    if ~iscell(new_kids)
      h = new_kids;
    else
      h  = [new_kids{:}];
    end
    
  otherwise
    warning('\nadd_errorbars: Unsupported object type: %s\n', obj_type)
    return
end

if size(errordata,2) ~= ncond
  error(sprintf('Number of conditions (columns) does not match: expected %d, found %d', ncond,size(errordata,2)))
end

if size(get(h(1),'XData'),2) ~= size(errordata,1)
  error(sprintf('Number of rows does not match: expected %d, found %d', size(get(h(1),'XData'),2), size(errordata,1)))
end

figure(gcf)
axes(gca)
hold on

for icond = 1:ncond
  xdata = get(h(icond),'XData');
  ydata = get(h(icond),'YData');
  
  nbar = size(xdata,2);
  
  for ibar = 1:nbar
    switch orientation
      case 'vertical'
        ylength = errordata(ibar,icond);
        cap_length = diff(xdata([2 3], ibar))/2;

        % Draw the vertical section
        xoffset = mean(xdata([2 3],ibar));
        ystart = ydata(2,ibar);
        if ystart >= 0
          ystop = ystart + ylength;
        else
          ystop = ystart - ylength;
        end
        if ~isempty(linecolor)
          line([xoffset xoffset], [ystart ystop],'color',linecolor);
        else
          line([xoffset xoffset], [ystart ystop]);
        end

        % Add the cap
        xstart = xoffset-cap_length/2;
        xstop = xstart + cap_length;
        if ~isempty(linecolor)
          line([xstart xstop], [ystop ystop],'color',linecolor);
        else
          line([xstart xstop], [ystop ystop]);
        end
      case 'horizontal'
        xlength = errordata(ibar,icond);
        cap_length = diff(ydata([2 3], ibar))/2;

        % Draw the horizontal section
        yoffset = mean(ydata([2 3],ibar));
        xstart = xdata(2,ibar);
        if xstart >= 0
          xstop = xstart + xlength;
        else
          xstop = xstart - xlength;
        end
        if ~isempty(linecolor)
          line([xstart xstop], [yoffset yoffset],'color',linecolor);
        else
          line([xstart xstop], [yoffset yoffset]);
        end

        % Add the cap
        ystart = yoffset-cap_length/2;
        ystop = ystart + cap_length;
        if ~isempty(linecolor)
          line([xstop xstop], [ystart ystop],'color',linecolor);
        else
          line([xstop xstop], [ystart ystop]);
        end

    end % switch
  end % for ibar=
end % for icond=

hold off
