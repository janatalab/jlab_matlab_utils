function sh = sph_elec(x,y,z,n);
%sh = sph_elec(x,y,z,n)
if size(x,1) ~= 1
	error('too many rows')
end;
	[az,el,r] = cart2sph(x,y,z);
	el = (pi/2)*ones(size(el,1),size(el,2)) - el;
	ass = legendre(n,cos(el));
	sh = zeros(2*n+1,size(el,2));
	for i=-n:n
		if i >= 0
			sh(i+n+1,:) = ass(abs(i)+1,:).*cos(abs(i)*(az));
		else
			sh(i+n+1,:) = ass(abs(i)+1,:).*sin(abs(i)*(az));
	end;
end;
