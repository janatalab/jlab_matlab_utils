function masked_map= mask_meg_fixval(x,y,z,map,theta_max,fixvalue)

[az,el,r] = cart2sph(x,y,z);
el = (pi/2)*ones(size(el,1),size(el,2)) - el;
masked_map = map;
for i=1:size(map,1)
	for j=1:size(map,2)
		if el(i,j) >= theta_max*pi/180;
			masked_map(i,j) = fixvalue;
		end
		if (abs(az(i,j)) < 0.7 & el(i,j) >= 80*pi/180)
			masked_map(i,j) = fixvalue;
		end;
	end;
end;


