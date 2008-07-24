function y = logspace(d1, d2, n)
%LOGSPACE Logarithmically spaced vector.
%   LOGSPACE(X1, X2) generates a row vector of 50 logarithmically
%   equally spaced points between decades 2^X1 and 2^X2.  If X2
%   is pi, then the points are between 2^X1 and pi.
%
%   LOGSPACE(X1, X2, N) generates N points.
%   For N < 2, LOGSPACE returns 2^X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   See also LINSPACE, :.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.11.4.1 $  $Date: 2004/07/05 17:01:21 $
%
%   5/31/07 S.T. adapted from logspace.m to fit a log2 curve

if nargin == 2
    n = 50;
end
n = double(n);

if d2 == pi
    d2 = log2(pi);
end

if d2 == single(pi)   
    d2 = log2(single(pi));
end

y = (2).^ [d1+(0:n-2)*(d2-d1)/(floor(n)-1), d2];
