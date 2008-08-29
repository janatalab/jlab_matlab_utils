function sh = sphharm(x,y,z,n,m);
%sh = sphharm(x,y,z,n,m);
[az,el,r] = cart2sph(x,y,z);
el = (pi/2)*ones(size(el,1),size(el,2)) - el;

for i = 1:size(el,1)
	ass = legendre(n,cos(el(i,:)));
	if m >= 0
	sh(i,:) = ass(m+1,:).*cos(m*(az(i,:)));
	else
	sh(i,:) = -ass(abs(m)+1,:).*sin(m*(az(i,:)));
	end;
end;
