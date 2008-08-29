function [xs,ys,zs] = polgrid(ns,ruu);
%[xs,ys,zs] = polgrid(ns,ruu);
%
%


xs = zeros(ns,ns);
ys = zeros(ns,ns);
zs = zeros(ns,ns);

for j = 1:ns
	xx = (j-(ns+1)/2)*240/(ns+1);
	for jj = 1:ns
		yy = (jj-(ns+1)/2)*240/(ns+1);
		th = sqrt(xx.^2+yy.^2)*pi/180;
		if xx == 0
			if yy >= 0 
				phi = pi/2;
			elseif yy < 0
				phi = -pi/2;
			else
				phi = 0;
			end;
		else
			phi = -atan(-yy/xx);
			if (yy <= 0 & xx <=0)
				phi = phi - pi;
			end;
			if (yy >= 0 & xx <= 0)
				phi = phi+pi;
			end;
			if (yy == 0 & xx >= 0)
				phi = 0;
			end;
			if (yy == 0 & xx <= 0)
				phi = pi;
			end;
		end;
		xs(j,jj) = ruu*sin(th)*cos(phi);
		ys(j,jj) = ruu*sin(th)*sin(phi);
		zs(j,jj) = ruu*cos(th);
	end;
end;
		
				


