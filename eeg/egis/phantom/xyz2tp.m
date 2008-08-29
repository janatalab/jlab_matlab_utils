function [tp,vtp] = xyz2tp(x,y,z,v)
%[tp,vtp] = xyz2tp(x,y,z,v)
tp = zeros(size(x,2)*size(x,1),2);
vtp = zeros(size(x,2)*size(x,1),1);
[az,el,r] = cart2sph(x,y,z);
el2 = pi/2*ones(size(el)) - el;
ic = 1;
for i = 1:size(az,1)
	for j = 1:size(az,2)
		tp(ic,1) = el2(i,j);
		tp(ic,2) = az(i,j);
		ic = ic+1;
		if nargin == 4
			vtp(ic) = v(i,j);
		end
	end
end
