function p = legendre_r(z,nmax);
p = zeros(nmax,size(z,2));
p(1,:) = z;
p(2,:) = 0.5*(3*z.^2 - ones(1,size(z,2)));
for n = 3:nmax
	p(n,:) = (z.*p(n-1,:)*(2*n-1)-(n-1)*p(n-2,:))/n;
end;

	
