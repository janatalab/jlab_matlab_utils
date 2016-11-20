function [x,y,idx,item] = spiral_layout(items)
% Generates a spiral layout on a square grid given a list of arbitrary length
% Generated are x and y coordinates, assuming that the starting coordinate
% is 0,0 with integer steps in positive and negative directions.
%
% Currently, only a counterclockwise spiral is generated, with the first
% step to the right from center


% 19Nov2016 Petr Janata

% Determine the number of elements
items = items(:);
nitems = length(items);

% Generate a square containing the next large square of nitems
dimlen = ceil(sqrt(nitems));

% Make sure that dimlen is odd
if ~mod(dimlen,2)
  dimlen = dimlen + 1;
end

npos = dimlen^2;

% Generate x and y axes with the starting position in the center
xpos = -(dimlen-1)/2:(dimlen-1)/2;
ypos = xpos;

max_x = 0;
max_y = 0;
min_x = 0;
min_y = 0;
curr_x = 0;
curr_y = 0;
xdir = 'inc';
ydir = 'same';
x = nan(npos,1);
y = nan(npos,1);

for ipos = 1:npos
  % Assign the current value
  x(ipos) = curr_x;
  y(ipos) = curr_y;
  
  % Update curr_x
  if strcmp(xdir,'inc')
    curr_x = curr_x + 1;
  elseif strcmp(xdir, 'dec')
    curr_x  = curr_x - 1;
  end
  
  % Update curr_y
  if strcmp(ydir, 'inc')
    curr_y = curr_y + 1;
  elseif strcmp(ydir, 'dec')
    curr_y = curr_y - 1;
  end
  
  % Now adjust the direction based on the updated curr_x and curr_y
  % This set of if-else statements generates a counterclockwise spiral
  if strcmp(xdir, 'inc') && curr_x > max_x
    xdir = 'same';
    ydir = 'inc';
    max_x = curr_x;
  elseif strcmp(ydir, 'inc') && curr_y > max_y
    ydir = 'same';
    xdir = 'dec';
    max_y = curr_y;
  elseif strcmp(xdir, 'dec') && curr_x < min_x
    xdir = 'same';
    ydir = 'dec';
    min_x = curr_x;
  elseif strcmp(ydir, 'dec') && curr_y < min_y
    xdir = 'inc';
    ydir = 'same';
    min_y = curr_y;
  end
  
end % for ipos=1:npos

idx = (1:npos)';

item = items;

if isnumeric(item)
  item(end+1:npos) = NaN;
elseif iscell(item)
  item(end+1:npos) = {''};
end

end % function
