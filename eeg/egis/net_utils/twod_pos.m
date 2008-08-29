function twod_dist = twod_pos(xelec,yelec,zelec,radius)

[az,el,r] = cart2sph(xelec,yelec,zelec);

el = pi*ones(size(el,1),size(el,2))/2 - el;

icount = 1
for i = 1:size(xelec,2)
	for j = i:size(xelec,2)
		zeroptaz = az(i);
		zeroptel = el(i);
		azgood = az(j);
		elgood = el(j);

		abgdist = acos(cos(elgood).*cos(zeroptel)+sin(elgood).*sin(zeroptel).*cos(azgood - zeroptaz));

		abgdist = abgdist*radius;

		twod_dist(icount) = abgdist;
		icount = icount+1;
	end;
end;


	
