function [wvf] = addramp(x,npts,type)
% [wvf] = addramp(x,npts,type);
%
% type:
%  'lin':  linear

% 03/09/00 PJ
% 02/07/05 PJ -- modified to handle matrices

msg = nargchk(2,3,nargin);
if msg
  error(msg)
end

if nargin < 3
  type = 'lin';
end

switch type
  case {'lin','linear'}
    ramp_vals = (0:1/(npts-1):1)';
  otherwise
    error(sprintf('Unknown ramp type: %s', type))
end

%x = x(:);

x(1:npts,:) = x(1:npts,:).*repmat(ramp_vals,1,size(x,2));  % onset ramp
x = flipud(x);
x(1:npts,:) = x(1:npts,:).*repmat(ramp_vals,1,size(x,2));  % offset ramp
wvf = flipud(x);
