function [k,kinv,a,ainv,e]= k_and_e(w,x,y,z)

n = size(x,2);
mm = 10;
for i=1:n
	if x(i) == 0
		x(i) = 0.001;
	end;
	if y(i) == 0
		y(i) = 0.001;
	end;
	if z(i) == 0
		z(i) = 0.001;
	end;
end;

for i=1:n
	for j=1:n
		s = x(i) - x(j);
		t = y(i) - y(j);
		r = z(i) - z(j);
		str = s.^2 + t.^2 + r.^2;
		k(i,j) = ((str+ w).^2)*log(str+w);
	end;
end;

for i=1:n
e(i,1) = 1;
e(i,2) = x(i);
e(i,3) = y(i);
e(i,4) = x(i).^2;
e(i,5) = x(i)*y(i);
e(i,6) = y(i).^2;
e(i,7) = z(i);
e(i,8) = z(i)*x(i);
e(i,9) = z(i)*y(i);
e(i,10) = z(i).^2;
end;
kinv = inv(k);

ke = kinv*e;

et = e';

a = et*ke;

ainv = inv(a);

