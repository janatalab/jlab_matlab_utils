function [a,b,c,d,e,f,g] = fender_4(nmax,s12,s13,r1,r2,r3,r4,rz)
% [a,b,c,d,e,f,g] = fender_4(nmax,s12,s13,r1,r2,r3,r4,rz)
% nmax = maximum n to calculate (a-g are sized to(1, nmax)
% s12 = conductivity brain/CSF
% s13 = conductivity brain/skull
% r1 = radii brain
% r2 = radii CSF
% r3 = radii skull
% r4 = radii scalp
% rz = radii source
a = zeros(1,nmax);
b = zeros(1,nmax);
c = zeros(1,nmax);
d = zeros(1,nmax);
e = zeros(1,nmax);
f = zeros(1,nmax);
g = zeros(1,nmax);

q = 1./4/pi/r1^2;

s34 = 1/s13;

s23 = s13/s12;

r12 = r1/r2;

r21 = r2/r1;

r23 = r2/r3;

r32 = r3/r2;

r43 = r4/r3;

r34 = r3/r4;

rz1 = rz/r1;

for n = 1:nmax
xn = n;
np = (n+1)/n;

vnfac = (r34^n - r43^(n+1))/(np*r34^n+r43^(n+1));
vn = (s34/np-vnfac)/(s34+vnfac);
wn = ((r23^n)/np - vn*r32^(n+1))/((r23^n)+vn*r32^(n+1));
yn = (s23/np-wn)/(s23+wn);
zn = (r12^n-np*yn*r21^(n+1))/(yn*r21^(n+1)+r12^n);
a(n) = q*((rz1^(n-1))*(n*zn+s12*(n+1)))/(s12-zn);
b(n) = (a(n)+xn*q*(rz1^(n-1)))/(yn*(r21^(n+1))+r12^n);
c(n) = yn*b(n);
d(n) = (b(n)+c(n))/(r23^n+vn*(r32^(n+1)));
e(n) = vn*d(n);
g(n) = (d(n)+e(n))/(np*r34^n+r43^(n+1));
f(n) = np*g(n);
end;
