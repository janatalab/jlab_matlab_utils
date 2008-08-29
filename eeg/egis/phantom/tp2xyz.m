function [x,y,z,v] = tp2xyz(tp,dim,vtp)
if dim == 1
	x = zeros(size(tp,1),1);
	y = zeros(size(tp,1),1);
	z = zeros(size(tp,1),1);
	v = zeros(size(tp,1),1);
	az = tp(:,2);
	el = pi/2*ones(size(tp,1),1) - tp(:,1);
	r = 9*ones(size(tp,1),1);
	[x,y,z] = sph2cart(az,el,r);
	if nargin == 3
	v = vtp;
		end;
elseif dim == 2
	x = zeros(sqrt(size(tp,1)),sqrt(size(tp,1)));
	y = zeros(sqrt(size(tp,1)),sqrt(size(tp,1)));
	z = zeros(sqrt(size(tp,1)),sqrt(size(tp,1)));
	v = zeros(sqrt(size(tp,1)),sqrt(size(tp,1)));
	az = tp(:,2);
	el = pi/2*ones(size(tp,1),1) - tp(:,1);
	r = 9*ones(size(tp,1),1);
	[xt,yt,zt] = sph2cart(az,el,r);
	ic = 1;
	for i = 1:size(x,1)
		for j = 1:size(x,2)
			x(i,j) = xt(ic);
			y(i,j) = yt(ic);
			z(i,j) = zt(ic);
			ic = ic +1;
			if nargin == 3
			v(i,j) = vtp(ic);
			end;
		end;
	end;
end;
