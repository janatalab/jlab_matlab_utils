function h = scatter_w_stats(xdata,ydata,params)
% Generates scatter plot and includes regression lines and statistics
%
% h = scatter_w_stats(xdata,ydata,params)

% 12Feb2015 Petr Janata

if ~all(size(xdata) == size(ydata))
  error('xdata and ydata must be the same size')
end

if size(xdata,2) == 1
  xdata = xdata(:);
  ydata = ydata(:);
end

ncols = size(xdata,2);

% Generate the scatter plot
for icol = 1:ncols
  h(icol) = scatter(xdata(:,icol), ydata(:,icol));
  hold 'on'
end

% Do the stats
for icol = 1:ncols
  [b,~,~,~,stats] = regress(ydata(:,icol),[ones(size(xdata(:,1))) xdata(:,icol)]);
  
  xlims = get(gca,'xlim');
  l(icol) = line(xlims,xlims*b(2)+b(1),'Color',get(h(icol),'CData'));
  
  Rsqr = stats(1);
  r(icol) = sqrt(Rsqr)*sign(b(2));
  p(icol) = stats(3);
  
  % Add the text
  if sign(b(2)) > 0
    horizontalAlign = 'right';
  else
    horizontalAlign = 'left';
  end
  verticalAlign = 'bottom';
  
  t(icol) = text(mean(xlims), mean(xlims)*b(2)+b(1), sprintf('r(%d) = %.2f, p = %.4f', ...
    sum(~isnan(ydata(:,icol)))-1, r(icol), p(icol)), ...
    'horizontalAlign', horizontalAlign, ...
    'verticalAlign', verticalAlign);
end

end