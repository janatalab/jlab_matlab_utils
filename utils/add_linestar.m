function add_linestar(cond_locs, cond_data, pval, params)

if length(cond_locs) > 2
  fprintf('Cannot plot lines between more than 2 conditions\n')
  return
end

if nargin < 4
  params = struct();
end

% Check to see if format for string has been specified
if ~isfield(params,'use_prob_text')
  params.use_prob_text = false;
end

if ~isfield(params, 'bottom_offset')
  params.bottom_offset = 1.1;
end

if ~isfield(params, 'top_offset')
  params.top_offset = 1.25;
end

% calculate where the lines will start at the bottom

bottom_points = cond_data*params.bottom_offset;
top_points = cond_data*params.top_offset;

stop_y = max(top_points);

% Plot vertical line segment above 1st condition
for icond = 1:2
  curr_cond = cond_locs(icond);
  start_x = curr_cond;
  stop_x = start_x;

  start_y = bottom_points(icond);

  line([start_x stop_x], [start_y stop_y], 'color', [0 0 0])
end

% Add the top horizontal line
line(cond_locs, [1 1]*stop_y, 'color',[0 0 0])

if params.use_prob_text
  str = sprintf('pval: %1.4f', pval);
else
  str = prob2str(pval,0.05,'*');
end
% Add significance as stars
text(mean(cond_locs), stop_y*1.00, str,  ...
  'horizontalalign', 'center', ...
  'verticalalign', 'bottom', ...
  'fontsize', 14)

