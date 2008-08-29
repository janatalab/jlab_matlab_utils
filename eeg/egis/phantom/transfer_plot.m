function [h,a,l] = transfer_plot(nmax,s12,s13,r1,r2,r3,r4,rz)
%[h,a,l] = transfer_plot(nmax,s12,s13,r1,r2,r3,r4,rz)

[a,b,c,d,e,f,g] = fender_4(nmax,s12,s13,r1,r2,r3,r4,rz);
h = f+g;

for i = 1:nmax
l(i) = i*(i+1)*h(i)/r4/r4;
end;
for i = 1:size(a,2)
a(i) = a(i)+i*((rz/r1)^(i-1))/(4*pi*r1^2);
end;
h = 4*pi*h./(2*[1:nmax]+ones(1,nmax));
a = 4*pi*a./(2*[1:nmax]+ones(1,nmax));
l = 4*pi*l./(2*[1:nmax]+ones(1,nmax));
