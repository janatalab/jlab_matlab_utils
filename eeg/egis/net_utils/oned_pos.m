function oned_array = oned_pos(good_chan,xelec,yelec,zelec,radius,calc_type)

[az,el,r] = cart2sph(xelec,yelec,zelec);

el = pi*ones(size(el,1),size(el,2))/2 - el;

if calc_type == 's'
		
	oned_array = az(good_chan) - az(good_chan(1))*ones(size(good_chan,1),size(good_chan,2));

elseif calc_type == 'g'
	
	zeroptaz = az(good_chan(1))*ones(size(good_chan,1),size(good_chan,2));
	zeroptel = el(good_chan(1))*ones(size(good_chan,1),size(good_chan,2));

	azgood = az(good_chan);
	elgood = el(good_chan);

	abgdist = acos(cos(elgood).*cos(zeroptel)+sin(elgood).*sin(zeroptel).*cos(azgood - zeroptaz));

	abgdist = abgdist*radius;

else

	error('unknown array type');

end;


	

 