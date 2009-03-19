function [histCenters,histRadii] = rose_scaled(varargin)
% plots an angle histogram such as the 'rose' function, but accepts scaled values
%
%   ROSE(THETA) plots the angle histogram for the angles in THETA.  
%   The angles in the vector THETA must be specified in radians.
%
%   See also ROSE, HIST, POLAR, COMPASS.
%
%   Clay M. Thompson 7-9-91
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.14.4.4 $  $Date: 2005/04/28 19:56:53 $
%
%
%   Modified 3/19/2009 Stefan Tomic
%     - Accepts different vector strength contributions to the histogram. 
%     - Added line style functionality and a flag whether or not toplot. 
%     - The output now is also a list of histogram angle centers
%       and radiuses, rather than the line data that the 'rose' function outputs.
%
%   New format is [histCenters,histRadii] = rose_scaled(angles,strengths,lineStyle,do_plot)
%     where
%       - strengths (optional) is a list of radii that correspond to
%         the angles list. Default strength = 1.
%       - lineStyle (optional) corresponds the the possible line styles given to
%         the plot function. Default is 'b-'
%       - do_plot (optional) is a true or false flag that dictates whether or
%         not to plot the histogram (default = true)


[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,3,nargs,'struct'));

theta = args{1};
if nargs > 1, 
  strengths = args{2}; 
else
  strengths = ones(size(theta));
end
if nargs > 2
  x = args{3};
end
if nargs > 3
  lineStyle = args{4};
else
  lineStyle = 'b-';
end

if nargs > 4
  do_plot = args{5};
else
  do_plot = 1;
end

if ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end
theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi
if (nargs <= 2),
  x = (0:19)*pi/10+pi/20;

elseif nargs==2,
  if ischar(strengths)
    error(id('NonNumericInput'),'Input arguments must be numeric.');
  end
  
elseif nargs==3
  if ischar(x)
    error(id('NonNumericInput'),'Input arguments must be numeric.');
  end

  if length(x)==1,
    x = (0:x-1)*2*pi/x + pi/x;
  else
    x = sort(rem(x(:)',2*pi));
  end

end
if ischar(x) || ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end

% Determine bin edges and get histogram
edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];
[dummy,bins] = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
%nn(end-1) = nn(end-1)+nn(end);
%nn(end) = [];

nBins = length(x);
for iBin = 1:nBins
  nn(iBin) = sum(strengths(bins == iBin));
end

% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);

if (do_plot)
  if ~isempty(cax)
    h = polar(cax,t,r,lineStyle);
  else
    h = polar(t,r,lineStyle);
  end
  
end

histCenters = x;
histRadii = nn;

function str=id(str)
str = ['MATLAB:rose:' str];
